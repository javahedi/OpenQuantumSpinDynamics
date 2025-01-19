# src/utils.jl
function setup_logging()
    logger = SimpleLogger(stdout, Logging.Info)
    global_logger(logger)
end

function summarize_results(results)
    success_count = count(results)
    failure_count = length(results) - success_count
    @info "Job summary: $success_count successful, $failure_count failed."
end