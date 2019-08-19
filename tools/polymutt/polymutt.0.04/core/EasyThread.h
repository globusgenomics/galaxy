#ifndef __EASYTHREAD_H__
#define __EASYTHREAD_H__

#define __THREAD_RUNNING__    0x1234
#define __THREAD_AVAILABLE__  0x8888

struct ThreadParameters
   {
   int  *  status;
   void ** result;
   void * (* func)(void *);
   void * arg;
   };

class ThreadManager
   {
   public:
      int   * threadStatus;
      void ** returnValues;
      ThreadParameters * threadParameters;
      int     maxThreads;

   ThreadManager(int threads);
   ~ThreadManager();

   int StartThread(void * (* func) (void *), void * arg);

   void WaitForAll();
   bool isIdle();

   private:
      int GetSlot();
#ifdef __WIN32__
      static void RunThread(void * runParameters);
#else
      static void * RunThread(void * runParameters);
#endif
   };

#endif


