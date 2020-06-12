#include "def.h"

uint64 NchoosK(uint32 n, uint32 k)
{
    if (k > n)
    {
        return 0;
    }
    uint64 r = 1;
    for (uint64 d = 1; d <= k; ++d)
    {
        r *= n--;
        r /= d;
    }
    return r;
}

struct JobData
{
public:
    uint32 jobId;        // multi threading
    uint64 testToDo;     // number of test to be performed in this job
    uint64 firstTest;    // initial state (processing start point)
    varIdx v[MAX_ORDER]; // initial state (processing start point)
    uint32 task;         // initial state (processing start point)
    uint32 taskOrder;    // initial state (processing start point)
    uint64 testCounter;  // used as global variable for recursive iteratation
    bool init;           // used as global variable for recursive iteratation
    bool done;           // used as global variable for recursive iteratation

    void Start()
    {
        testCounter = 0;
        init = true;
        done = false;
    }
    void Next()
    {
        init = false;
        testCounter++;
        if (testCounter == testToDo)
        {
            done = true;
        }
    }
    void Print()
    {
        printf("\n%15llu%15llu%8u%4u [", firstTest, testToDo, task, taskOrder);
        for (uint32 i = 0; i < taskOrder; i++)
            printf("%8u", v[i]);
        printf("]");
    }
};

class WorkLoad
{
public:
    uint32 numSets;
    uint32 *setsLen;
    varIdx **sets;

    uint32 numTasks;
    uint32 *tasksLen;
    uint32 **tasks;

    uint64 *numTestPerTask;
    uint64 numTest;

    uint32 numJobs;
    uint32 firstJob;
    uint32 lastJob;
    uint32 jobsToDo;

    JobData *jobData;
    uint32 jobCounter;  // used as global variable for recursive iteratation
    uint64 testCounter; // used as global variable for recursive iteratation

    // rnl: Recursive Nested Loop
    void RnlWorkLoadDivider(uint32 task, varIdx *v, uint32 level)
    {
        uint32 *taskStes = tasks[task];
        uint32 taskOrder = tasksLen[task];

        uint start;
        if ((level > 0) && (taskStes[level] == taskStes[level - 1]))
        {
            start = v[level - 1] + 1;
        }
        else // ((level > 0) || (taskStes[level] != taskStes[level - 1]))
        {
            start = 0;
        }

        for (v[level] = start; v[level] < setsLen[taskStes[level]]; v[level]++)
        {
            if (taskOrder > (level + 1))
            {
                RnlWorkLoadDivider(task, v, level + 1);
            }
            else // if last level
            {
                if (testCounter == jobData[jobCounter].firstTest)
                {
                    memcpy(jobData[jobCounter].v, v, sizeof(varIdx) * MAX_ORDER);
                    jobData[jobCounter].task = task;
                    jobData[jobCounter].taskOrder = taskOrder;
                    jobCounter++;
                }
                testCounter++;
            }
        }
    }

    void RnlPrintTestPerTask(uint32 task, varIdx *v, uint32 level)
    {
        uint32 *taskStes = tasks[task];
        uint32 taskOrder = tasksLen[task];

        uint start;
        if ((level > 0) && (taskStes[level] == taskStes[level - 1]))
        {
            start = v[level - 1] + 1;
        }
        else // ((level > 0) || (taskStes[level] != taskStes[level - 1]))
        {
            start = 0;
        }

        for (v[level] = start; v[level] < setsLen[taskStes[level]]; v[level]++)
        {
            if (taskOrder > (level + 1))
            {
                RnlPrintTestPerTask(task, v, level + 1);
            }
            else // if last level
            {
                printf("\n");
                for (uint32 i = 0; i < taskOrder; i++)
                    printf("%4u%c", v[i], 'A' + taskStes[i]);
            }
        }
        return;
    }

    void RnlPrintTestPerJob(uint32 task, uint32 job, uint32 level)
    {
        uint32 *taskStes = tasks[task];
        uint32 taskOrder = tasksLen[task];
        varIdx *v = jobData[job].v;

        uint start;
        if (jobData[job].init)
        {
            start = v[level];
        }
        else if ((level > 0) && (taskStes[level] == taskStes[level - 1]))
        {
            start = v[level - 1] + 1;
        }
        else // ((level > 0) || (taskStes[level] != taskStes[level - 1]))
        {
            start = 0;
        }

        for (v[level] = start; v[level] < setsLen[taskStes[level]]; v[level]++)
        {
            if (jobData[job].done)
                break;

            if (taskOrder > (level + 1))
            {
                RnlPrintTestPerJob(task, job, level + 1);
            }
            else // if last level
            {
                printf("\n task: %8u", task);
                for (uint32 i = 0; i < taskOrder; i++)
                    printf("%4u%c", v[i], 'A' + taskStes[i]);

                jobData[job].Next();
                if (jobData[job].done)
                    break;
            }
        }
        return;
    }

    void RnlPerformJob(uint32 task, uint32 job, uint32 level)
    {
        uint32 *taskStes = tasks[task];
        uint32 taskOrder = tasksLen[task];
        varIdx *v = jobData[job].v;

        uint start;
        if (jobData[job].init)
        {
            start = v[level];
        }
        else if ((level > 0) && (taskStes[level] == taskStes[level - 1]))
        {
            start = v[level - 1] + 1;
        }
        else // ((level > 0) || (taskStes[level] != taskStes[level - 1]))
        {
            start = 0;
        }

        for (v[level] = start; v[level] < setsLen[taskStes[level]]; v[level]++)
        {
            if (jobData[job].done)
                break;

            if (taskOrder > (level + 1))
            {
                // **** OR

                RnlPrintTestPerJob(task, job, level + 1);
            }
            else // if last level
            {
                //*** OR  CONTINGENCTY

                jobData[job].Next();
                if (jobData[job].done)
                    break;
            }
        }
        return;
    }

    void WorkLoadDivider()
    {
        jobCounter = 0;
        testCounter = 0;
        for (uint32 i = 0; i < numTasks; i++)
        {
            varIdx v[MAX_ORDER];
            memset(v, 0, sizeof(varIdx) * MAX_ORDER);
            RnlWorkLoadDivider(i, v, 0);
        }
        for (uint32 i = 0; i < numJobs; i++)
            jobData[i].Print();
    }

    void PrintTestPerTask()
    {
        jobCounter = 0;
        testCounter = 0;
        for (uint32 i = 0; i < numTasks; i++)
        {
            varIdx v[MAX_ORDER];
            memset(v, 0, sizeof(varIdx) * MAX_ORDER);
            printf("\n>>> Tests in task %u:", i);
            RnlPrintTestPerTask(i, v, 0);
        }
    }

    void PrintTestPerJob()
    {
        for (uint32 i = 0; i < numJobs; i++)
        {
            uint32 task = jobData[i].task;
            jobData[i].Start();

            printf("\n>>> Tests in job %u:", i);
            while (true)
            {
                RnlPrintTestPerJob(task, i, 0);
                if (jobData[i].done)
                    break;
                task++;
            }
        }
    }

    void Test()
    {
        numJobs = 4;
        firstJob = 2;
        lastJob = 3;
        jobsToDo = (lastJob - firstJob) + 1;
        numSets = 3;
        numTasks = 3;

        jobData = new JobData[numJobs];
        setsLen = new uint32[numSets];
        sets = new varIdx *[numSets];
        tasksLen = new uint32[numTasks];
        tasks = new uint32 *[numTasks];
        numTestPerTask = new uint64[numTasks];

        NULL_CHECK(jobData);
        NULL_CHECK(setsLen);
        NULL_CHECK(sets);
        NULL_CHECK(tasksLen);
        NULL_CHECK(tasks);
        NULL_CHECK(numTestPerTask);

        memset(jobData, 0, sizeof(JobData) * numJobs);

        setsLen[0] = 100;
        setsLen[1] = 3;
        setsLen[2] = 5;

        tasksLen[0] = 2;
        tasksLen[1] = 3;
        tasksLen[2] = 4;

        for (uint32 i = 0; i < numSets; i++)
        {
            sets[i] = new varIdx[setsLen[i]];
            NULL_CHECK(sets[i]);
        }

        for (uint32 i = 0; i < numTasks; i++)
        {
            tasks[i] = new varIdx[tasksLen[i]];
            NULL_CHECK(tasks[i]);
        }

        tasks[0][0] = 1;
        tasks[0][1] = 2;
        tasks[1][0] = 1;
        tasks[1][1] = 1;
        tasks[1][2] = 2;
        tasks[2][0] = 2;
        tasks[2][1] = 1;
        tasks[2][2] = 1;
        tasks[2][3] = 1;

        ComputeNumTest();
        printf("\n>>>%llu<<<", numTest);

        uint64 testPerJob = numTest / numJobs;
        uint64 testLastJob = testPerJob + (numTest % numJobs);

        for (uint32 i = 0; i < numJobs; i++)
        {
            jobData[i].firstTest = i * testPerJob;
            jobData[i].testToDo = testPerJob;
            if (i == (numJobs - 1))
                jobData[i].testToDo = testLastJob;
        }

        PrintTestPerTask();
        WorkLoadDivider();
        PrintTestPerJob();
    }

    ~WorkLoad()
    {

        for (uint32 i = 0; i < numSets; i++)
        {
            delete[] sets[i];
        }

        for (uint32 i = 0; i < numTasks; i++)
        {
            delete[] tasks[i];
        }

        delete[] jobData;
        delete[] setsLen;
        delete[] sets;
        delete[] tasksLen;
        delete[] tasks;
        delete[] numTestPerTask;
    }

    void ComputeNumTest()
    {
        numTest = 0;
        for (uint32 i = 0; i < numTasks; i++)
        {
            numTestPerTask[i] = 1;
            for (uint32 j = 0; j < tasksLen[i]; j++)
            {
                uint32 set = tasks[i][j];
                uint32 repeat = 1;
                while (true)
                {
                    if (j < (tasksLen[i] - 1))
                    {
                        if (tasks[i][j + 1] == set)
                        {
                            j++;
                            repeat++;
                        }
                        else
                        {
                            break;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
                numTestPerTask[i] *= NchoosK(setsLen[set], repeat);
            }
            numTest += numTestPerTask[i];
        }
        return;
    }

    void RunJob(uint32 job)
    {
        job--; // 1-base to 0-base index
        uint32 task = jobData[job].task;
        jobData[job].Start();

        printf("\n>>> Tests in job %u:", job);
        while (true)
        {
            RnlPrintTestPerJob(task, job, 0);
            if (jobData[job].done)
                break;
            task++;
        }
    }
};

struct ThreadData
{
    uint32 jobId;
    WorkLoad *workLoad;
};

void *threadFunction(void *data)
{
    ThreadData *td = (ThreadData *)data;
    td->workLoad->RunJob(td->jobId);
    return NULL;
}

void RunAllJobs(WorkLoad &workLoad)
{

    ThreadData *tds = new ThreadData[workLoad.jobsToDo];
    NULL_CHECK(tds);

    pthread_t *threads = new pthread_t[workLoad.jobsToDo];
    NULL_CHECK(threads);

    for (uint32 t = 0; t < workLoad.jobsToDo; t++)
    {
        tds[t].jobId = t + workLoad.firstJob;
        tds[t].workLoad = &workLoad;
        pthread_create(&threads[t], NULL, threadFunction, &tds[t]);
        pthread_join(threads[t], NULL);
    }
    // for (uint32 t = 0; t < workLoad.jobsToDo; t++)
    // {
    //     pthread_join(threads[t], NULL);
    // }

    delete[] threads;
    delete[] tds;
}