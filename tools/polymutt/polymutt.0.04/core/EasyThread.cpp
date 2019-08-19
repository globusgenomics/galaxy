#include "EasyThread.h"

#include <stdlib.h>

#ifdef __WIN32__
#include <process.h>
#else
#include <pthread.h>
#endif

ThreadManager::ThreadManager(int threads)
   {
   maxThreads = threads > 1 ? threads : 1;

   returnValues = new void * [maxThreads];
   threadStatus = new int [maxThreads];
   threadParameters = new ThreadParameters [maxThreads];

   for (int i = 0; i < maxThreads; i++)
      {
      threadStatus[i] = __THREAD_AVAILABLE__;
      returnValues[i] = NULL;
      }
   }

ThreadManager::~ThreadManager()
   {
   delete [] threadParameters;
   delete [] threadStatus;
   delete [] returnValues;
   }

int ThreadManager::StartThread(void * (* func) (void *), void * arg)
   {
   if (maxThreads > 1)
      {
      int slot = GetSlot();

      threadStatus[slot] = __THREAD_RUNNING__;

      threadParameters[slot].func = func;
      threadParameters[slot].arg  = arg;
      threadParameters[slot].status = &threadStatus[slot];
      threadParameters[slot].result = &returnValues[slot];

#ifdef __WIN32__
      while (_beginthread(RunThread, 32 * 1024, (void *)(threadParameters + slot))
             == (unsigned long)-1) ;
#else
      RunThread((void *) threadParameters);
#endif

      return slot;
      }
   func(arg);
   return 1;
   }

#ifdef __WIN32__
void ThreadManager::RunThread(void * arg)
#else
void * ThreadManager::RunThread(void * arg)
#endif
   {
   ThreadParameters * tp = (ThreadParameters *) arg;

   *(tp->result) = tp->func(tp->arg);

   *(tp->status) = __THREAD_AVAILABLE__;

#ifndef __WIN32__
   return NULL;
#endif
   }

int ThreadManager::GetSlot()
   {
   while ( true )
      {
      for (int i = 0; i < maxThreads; i++)
         if (threadStatus[i] == __THREAD_AVAILABLE__)
            return i;

      // waste some time ?
      ; ; ;
      }
   };

bool ThreadManager::isIdle()
   {
   for (int i = 0; i < maxThreads; i++)
      if (threadStatus[i] == __THREAD_RUNNING__)
         return false;
   return true;
   }

void ThreadManager::WaitForAll()
   {
   while (!isIdle())
      {
      // waste some time ?
      ; ; ;
      }
   }

