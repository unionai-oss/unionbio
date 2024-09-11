from flytekit import task, workflow, LaunchPlan, CronSchedule

@task
def get_batch() -> List[Jobs]:
    # return n oldest jobs from db in "TODO" state
    return query_jobs()

@task
def run_af(jobs: List[Jobs]) -> List[Results]:
    results = []
    for job in jobs:
        result = fold(job)
        results.append(result)
        update_job_in_db(job, result)
    return results

@workflow
def af_wf() -> List[Results]:
    jobs = get_batch()
    return run_af(jobs)

af_lp = LaunchPlan.create("af_wf", af_wf, schedule=CronSchedule(schedule="0 * * * *"))

