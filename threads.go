package gravitree

import (
	"fmt"
	"runtime"
)

// WorkerQueue is (you guessed it) a simple worker queue
// implementation. 
//
// How to use:
//
// workers, jobs := 16, 1000
// WorkerQueue(workers, jobs, func(worker, job int) {
//     /* use resources associated with [worker] to do [job] */
// })
func WorkerQueue(workers, jobs int, work func(worker, job int)) {
	if workers > jobs { workers = jobs }
	
	jobChan := make(chan int, jobs)
	lockChan := make(chan int, workers)

	for i := 0; i < jobs; i++ { jobChan <- i }
	
	for i := 0; i < workers; i++ {
		go func(workerIdx int) {
			for j := range jobChan {
				work(workerIdx, j)
			}
			lockChan <- 0
		}(i)
	}

	close(jobChan)
	for i := 0; i < workers; i++ { <-lockChan }
}

var nWorkers = runtime.GOMAXPROCS(-1)

// SetThreads sets the number of threads to use for parallel
// calculations. If you set n to -1, one thread will be used
// for each CPU.
func SetThreads(n int) {
	if n == -1 {
		runtime.GOMAXPROCS(runtime.NumCPU())
	} else if n > 0 {
		runtime.GOMAXPROCS(n)
	} else {
		panic(fmt.Sprintf("Invalid thread count: %d", n))
	}

	nWorkers = runtime.GOMAXPROCS(-1)
}
