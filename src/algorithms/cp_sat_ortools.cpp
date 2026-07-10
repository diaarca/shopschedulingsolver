#include "shopschedulingsolver/algorithms/cp_sat_ortools.hpp"

#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/cp_model_solver.h"

#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

using namespace operations_research;
using namespace operations_research::sat;

namespace
{

struct OperationType
{
    IntVar start;
    IntVar end;
    IntervalVar interval;
};

}

shopschedulingsolver::Output shopschedulingsolver::cp_sat_ortools(
        const Instance& instance,
        const Parameters& parameters)
{
    Output output(instance);

    if (instance.objective() == Objective::Throughput) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "invalid 'objective'; "
                "objective Throughput is not supported.");
    }
    if (instance.flexible()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "invalid 'flexible'; "
                "flexible instances are not supported.");
    }
    if (instance.operations_arbitrary_order()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "invalid 'operations_arbitrary_order'; "
                "open-shop instances are not supported.");
    }
    if (instance.blocking()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "invalid 'blocking'; "
                "blocking instances are not supported.");
    }
    if (instance.no_idle()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "invalid 'no_idle'; "
                "no-idle instances are not supported.");
    }
    if (instance.mixed_no_idle()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "invalid 'mixed_no_idle'; "
                "mixed no-idle instances are not supported.");
    }
    if (instance.permutation()) {
        throw std::invalid_argument(
                FUNC_SIGNATURE + ": "
                "invalid 'permutation'; "
                "permutation instances are not supported.");
    }

    AlgorithmFormatter algorithm_formatter(instance, parameters, output);
    algorithm_formatter.start("CP-SAT (OR-Tools)");

    // Compute time horizon for variable bounds
    int64_t max_release_date = 0;
    int64_t total_duration = 0;
    for (JobId job_id = 0; job_id < instance.number_of_jobs(); ++job_id) {
        const auto& job = instance.job(job_id);
        max_release_date = std::max(max_release_date, job.release_date);
        for (const auto& operation: job.operations) {
            Time max_p = 0;
            for (const auto& alternative: operation.alternatives)
                max_p = std::max(max_p, alternative.processing_time);
            total_duration += max_p;
        }
    }
    int64_t horizon = max_release_date + total_duration;

    CpModelBuilder cp_model;

    std::vector<std::vector<OperationType>> all_operations(
        instance.number_of_jobs());
    std::vector<std::vector<IntervalVar>> machine_to_intervals(
        instance.number_of_machines());
    std::vector<IntVar> job_ends(instance.number_of_jobs());

    // Model variables (jobs and operations)
    for (JobId job_id = 0; job_id < instance.number_of_jobs(); ++job_id) {
        const auto& job = instance.job(job_id);
        all_operations[job_id].resize(job.operations.size());

        job_ends[job_id] = cp_model.NewIntVar({0, horizon})
                               .WithName("job_end_" + std::to_string(job_id));

        // Per operation variables (start, end, interval)
        for (OperationId operation_id = 0;
                operation_id < (OperationId)job.operations.size();
                ++operation_id) {
            const auto& operation = job.operations[operation_id];

            // We don't handle the flexible variant for now so we only care on
            // the only first machine alternative
            const auto& alternative = operation.alternatives[0];
            MachineId machine = alternative.machine_id;
            Time duration = alternative.processing_time;

            std::string suffix =
                "_" + std::to_string(job_id) + "_" + std::to_string(operation_id);
            IntVar start = cp_model.NewIntVar({job.release_date, horizon})
                               .WithName(std::string("start") + suffix);
            IntVar end = cp_model.NewIntVar({job.release_date, horizon})
                             .WithName(std::string("end") + suffix);
            IntervalVar interval =
                cp_model.NewIntervalVar(start, duration, end)
                    .WithName(std::string("interval") + suffix);

            all_operations[job_id][operation_id] = OperationType{start, end, interval};
            machine_to_intervals[machine].push_back(interval);
        }
    }

    // No overlap on each machine
    for (MachineId machine_id = 0;
            machine_id < instance.number_of_machines();
            ++machine_id) {
        cp_model.AddNoOverlap(machine_to_intervals[machine_id]);
    }

    // Job constraints
    for (JobId job_id = 0; job_id < instance.number_of_jobs(); ++job_id) {
        const auto& job = instance.job(job_id);
        std::vector<IntVar> op_ends;

        for (OperationId operation_id = 0;
                operation_id < (OperationId)job.operations.size();
                ++operation_id) {
            op_ends.push_back(all_operations[job_id][operation_id].end);
        }

        cp_model.AddMaxEquality(job_ends[job_id], op_ends);

        // Job's operations precedence constraints
        for (OperationId operation_id = 0;
                operation_id < (OperationId)job.operations.size() - 1;
                ++operation_id) {
            if (instance.no_wait()) {
                // No time between two consecutive job's operations
                cp_model.AddEquality(
                        all_operations[job_id][operation_id + 1].start,
                        all_operations[job_id][operation_id].end);
            } else {
                cp_model.AddGreaterOrEqual(
                        all_operations[job_id][operation_id + 1].start,
                        all_operations[job_id][operation_id].end);
            }
        }
    }

    // Objective function
    if (instance.objective() == Objective::Makespan) {
        IntVar makespan = cp_model.NewIntVar({0, horizon}).WithName("Makespan");
        std::vector<IntVar> active_job_ends;
        for (JobId job_id = 0; job_id < instance.number_of_jobs(); ++job_id)
            active_job_ends.push_back(job_ends[job_id]);
        cp_model.AddMaxEquality(makespan, active_job_ends);
        cp_model.Minimize(makespan);

    } else if (instance.objective() == Objective::TotalFlowTime) {
        LinearExpr total_flow_time;
        for (JobId job_id = 0; job_id < instance.number_of_jobs(); ++job_id) {
            const auto& job = instance.job(job_id);
            total_flow_time += job.weight * (job_ends[job_id] - job.release_date);
        }
        cp_model.Minimize(total_flow_time);

    } else if (instance.objective() == Objective::TotalTardiness) {
        LinearExpr total_tardiness;
        for (JobId job_id = 0; job_id < instance.number_of_jobs(); ++job_id) {
            const auto& job = instance.job(job_id);
            // Ignore the jobs that don't have deadlines
            if (job.due_date == -1)
                continue;
            IntVar tardiness = cp_model.NewIntVar({0, horizon})
                .WithName("tardiness_" + std::to_string(job_id));
            cp_model.AddGreaterOrEqual(tardiness, job_ends[job_id] - job.due_date);
            total_tardiness += job.weight * tardiness;
        }
        cp_model.Minimize(total_tardiness);
    }

    operations_research::sat::SatParameters cp_parameters;

    // Write solver output to file.
    // TODO
    //cp_parameters.set_log_search_progress(true);

    // Set time limit.
    cp_parameters.set_max_time_in_seconds(parameters.timer.remaining_time());

    // Solve the model
    operations_research::sat::CpSolverResponse response = operations_research::sat::SolveWithParameters(cp_model.Build(), cp_parameters);

    // Build the Output object
    if (response.status() == CpSolverStatus::OPTIMAL
            || response.status() == CpSolverStatus::FEASIBLE) {
        SolutionBuilder solution_builder;
        solution_builder.set_instance(instance);

        for (JobId job_id = 0; job_id < instance.number_of_jobs(); ++job_id) {
            const auto& job = instance.job(job_id);
            for (OperationId operation_id = 0;
                    operation_id < (OperationId)job.operations.size();
                    ++operation_id) {
                Time start = SolutionIntegerValue(
                    response, all_operations[job_id][operation_id].start);
                solution_builder.append_operation(job_id, operation_id, 0, start);
            }
        }
        solution_builder.sort_machines();
        solution_builder.sort_jobs();
        algorithm_formatter.update_solution(solution_builder.build(), "");

        switch (instance.objective()) {
        case Objective::Makespan:
            algorithm_formatter.update_makespan_bound(
                std::round(response.best_objective_bound()), "");
            break;
        case Objective::TotalFlowTime:
            algorithm_formatter.update_total_flow_time_bound(
                std::round(response.best_objective_bound()), "");
            break;
        case Objective::TotalTardiness:
            algorithm_formatter.update_total_tardiness_bound(
                std::round(response.best_objective_bound()), "");
            break;
        default:
            // shouldn't happens
            break;
        }
    }

    algorithm_formatter.end();
    return output;
}
