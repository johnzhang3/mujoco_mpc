// Copyright 2023 DeepMind Technologies Limited
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "grpc/agent_service.h"

#include <cstring>
#include <memory>
#include <string_view>
#include <vector>

#include <absl/log/check.h>
#include <absl/strings/str_format.h>
#include <grpcpp/server_context.h>
#include <grpcpp/support/status.h>
#include <mujoco/mujoco.h>
#include "grpc/agent.pb.h"
#include "grpc/grpc_agent_util.h"
#include "mjpc/task.h"

namespace agent_grpc {

using ::agent::GetActionRequest;
using ::agent::GetActionResponse;
using ::agent::GetStateRequest;
using ::agent::GetStateResponse;
using ::agent::InitRequest;
using ::agent::InitResponse;
using ::agent::PlannerStepRequest;
using ::agent::PlannerStepResponse;
using ::agent::ResetRequest;
using ::agent::ResetResponse;
using ::agent::SetCostWeightsRequest;
using ::agent::SetCostWeightsResponse;
using ::agent::SetStateRequest;
using ::agent::SetStateResponse;
using ::agent::SetTaskParametersRequest;
using ::agent::SetTaskParametersResponse;
using ::agent::StepRequest;
using ::agent::StepResponse;

mjpc::Task* task = nullptr;
// model used for physics
mjModel* model = nullptr;
// model used for planning, owned by the Agent instance.
mjModel* agent_model = nullptr;

void residual_sensor_callback(const mjModel* m, mjData* d, int stage) {
  // with the `m == model` guard in place, no need to clear the callback.
  if (m == agent_model || m == model) {
    if (stage == mjSTAGE_ACC) {
      task->Residual(m, d, d->sensordata);
    }
  }
}

mjModel* LoadModelFromString(std::string_view xml, char* error,
                             int error_size) {
  static constexpr char file[] = "temporary-filename.xml";
  // mjVFS structs need to be allocated on the heap, because it's ~2MB
  auto vfs = std::make_unique<mjVFS>();
  mj_defaultVFS(vfs.get());
  mj_makeEmptyFileVFS(vfs.get(), file, xml.size());
  int file_idx = mj_findFileVFS(vfs.get(), file);
  memcpy(vfs->filedata[file_idx], xml.data(), xml.size());
  mjModel* m = mj_loadXML(file, vfs.get(), error, error_size);
  mj_deleteFileVFS(vfs.get(), file);
  return m;
}

mjModel* LoadModelFromBytes(std::string_view mjb) {
  static constexpr char file[] = "temporary-filename.mjb";
  // mjVFS structs need to be allocated on the heap, because it's ~2MB
  auto vfs = std::make_unique<mjVFS>();
  mj_defaultVFS(vfs.get());
  mj_makeEmptyFileVFS(vfs.get(), file, mjb.size());
  int file_idx = mj_findFileVFS(vfs.get(), file);
  memcpy(vfs->filedata[file_idx], mjb.data(), mjb.size());
  mjModel* m = mj_loadModel(file, vfs.get());
  mj_deleteFileVFS(vfs.get(), file);
  return m;
}

grpc::Status AgentService::Init(grpc::ServerContext* context,
                                const InitRequest* request,
                                InitResponse* response) {
  std::string_view task_id = request->task_id();
  agent_.SetTaskList(tasks_);
  int task_index = agent_.GetTaskIdByName(task_id);
  if (task_index == -1) {
    return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT,
                        absl::StrFormat("Invalid task_id: '%s'", task_id));
  }

  agent_.SetTaskByIndex(task_index);
  // TODO(khartikainen): is this needed?
  agent_.gui_task_id = task_index;

  mjModel* tmp_model;
  char load_error[1024] = "";

  if (request->has_model() && request->model().has_mjb()) {
    std::string_view model_mjb_bytes = request->model().mjb();
    // TODO(khartikainen): Add error handling for mjb loading.
    tmp_model = LoadModelFromBytes(model_mjb_bytes);
  } else if (request->has_model() && request->model().has_xml()) {
    std::string_view model_xml = request->model().xml();
    tmp_model = LoadModelFromString(model_xml, load_error, sizeof(load_error));
  } else {
    tmp_model = mj_loadXML(agent_.ActiveTask()->XmlPath().c_str(), nullptr,
                           load_error, sizeof(load_error));
  }

  if (!tmp_model) {
    return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT,
                        absl::StrFormat("Model load error: '%s'", load_error));
  }

  agent_.Initialize(tmp_model);
  mj_deleteModel(tmp_model);
  agent_.Allocate();
  agent_.Reset();

  task = agent_.ActiveTask();
  CHECK_EQ(agent_model, nullptr)
      << "Multiple instances of AgentService detected.";
  agent_model = agent_.GetModel();
  // copy the model before agent model's timestep and integrator are updated
  CHECK_EQ(model, nullptr)
      << "Multiple instances of AgentService detected.";
  model = mj_copyModel(nullptr, agent_model);
  data_ = mj_makeData(model);
  mjcb_sensor = residual_sensor_callback;

  agent_.SetState(data_);

  agent_.plan_enabled = true;
  agent_.action_enabled = true;

  return grpc::Status::OK;
}

AgentService::~AgentService() {
  if (data_) mj_deleteData(data_);
  if (model) mj_deleteModel(model);
  model = nullptr;
  // no need to delete agent_model and task, since they're owned by agent_.
  agent_model = nullptr;
  task = nullptr;
  mjcb_sensor = nullptr;
}

grpc::Status AgentService::GetState(grpc::ServerContext* context,
                                    const GetStateRequest* request,
                                    GetStateResponse* response) {
  if (!Initialized()) {
    return {grpc::StatusCode::FAILED_PRECONDITION, "Init not called."};
  }
  agent_.ActiveState().CopyTo(model, data_);
  return grpc_agent_util::GetState(model, data_, response);
}

grpc::Status AgentService::SetState(grpc::ServerContext* context,
                                    const SetStateRequest* request,
                                    SetStateResponse* response) {
  if (!Initialized()) {
    return {grpc::StatusCode::FAILED_PRECONDITION, "Init not called."};
  }
  grpc::Status status =
      grpc_agent_util::SetState(request, &agent_, model, data_);
  if (!status.ok()) return status;

  mj_forward(model, data_);
  task->Transition(model, data_);

  return grpc::Status::OK;
}

grpc::Status AgentService::GetAction(grpc::ServerContext* context,
                                     const GetActionRequest* request,
                                     GetActionResponse* response) {
  if (!Initialized()) {
    return {grpc::StatusCode::FAILED_PRECONDITION, "Init not called."};
  }
  return grpc_agent_util::GetAction(request, &agent_, response);
}

grpc::Status AgentService::PlannerStep(grpc::ServerContext* context,
                                       const PlannerStepRequest* request,
                                       PlannerStepResponse* response) {
  if (!Initialized()) {
    return {grpc::StatusCode::FAILED_PRECONDITION, "Init not called."};
  }
  agent_.plan_enabled = true;
  agent_.PlanIteration(&thread_pool_);

  return grpc::Status::OK;
}

grpc::Status AgentService::Step(grpc::ServerContext* context,
                                const StepRequest* request,
                                StepResponse* response) {
  if (!Initialized()) {
    return {grpc::StatusCode::FAILED_PRECONDITION, "Init not called."};
  }
  mjpc::State& state = agent_.ActiveState();
  state.CopyTo(model, data_);
  if (request->use_previous_policy()) {
    return {grpc::StatusCode::UNIMPLEMENTED,
            "use_previous_policy not implemented."};
  }
  agent_.ActivePlanner().ActionFromPolicy(data_->ctrl, state.state().data(),
                                          state.time());
  mj_step(model, data_);
  state.Set(model, data_);
  return grpc::Status::OK;
}

grpc::Status AgentService::Reset(grpc::ServerContext* context,
                                 const ResetRequest* request,
                                 ResetResponse* response) {
  if (!Initialized()) {
    return {grpc::StatusCode::FAILED_PRECONDITION, "Init not called."};
  }
  return grpc_agent_util::Reset(&agent_, agent_.GetModel(), data_);
}

grpc::Status AgentService::SetTaskParameters(
    grpc::ServerContext* context, const SetTaskParametersRequest* request,
    SetTaskParametersResponse* response) {
  if (!Initialized()) {
    return {grpc::StatusCode::FAILED_PRECONDITION, "Init not called."};
  }
  return grpc_agent_util::SetTaskParameters(request, &agent_);
}

grpc::Status AgentService::SetCostWeights(
    grpc::ServerContext* context, const SetCostWeightsRequest* request,
    SetCostWeightsResponse* response) {
  if (!Initialized()) {
    return {grpc::StatusCode::FAILED_PRECONDITION, "Init not called."};
  }
  return grpc_agent_util::SetCostWeights(request, &agent_);
}

}  // namespace agent_grpc
