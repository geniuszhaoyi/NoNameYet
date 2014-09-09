/** @file
@brief Macros for compiling in both linux and windows.
@author Yi Zhao
@note Document are gathered in PREDEFINED __linux
*/
#ifdef  _WIN32
#endif // _WIN32
#ifdef  __linux
#include <pthread.h>
#include <semaphore.h>
#endif // __linux

#ifdef  _WIN32
typedef int mos_pthread_mutex_t;
typedef int mos_sem_t;
typedef int mos_pthread_t;
typedef int mos_sem_t;
typedef int pthread_attr_t;
#define mos_pthread_mutex_init(a,b) a;b
#define mos_sem_init(a,b,c) a;b;c
#define mos_pthread_mutex_lock(a) a
#define mos_pthread_mutex_unlock(a) a
#define mos_sem_wait(a) a
#define mos_sem_post(a) a
#define mos_pthread_join(a,b) a;b
#define mos_pthread_detach(a) a
#define mos_pthread_self() 0
#define mos_pthread_create(a,b,c,d) 0;*(a)=0;c(NULL)
#endif // __linux
#ifdef  __linux
typedef pthread_mutex_t mos_pthread_mutex_t;
typedef sem_t mos_sem_t;
typedef pthread_t mos_pthread_t;
typedef sem_t mos_sem_t;
typedef pthread_attr_t mos_pthread_attr_t;
#define mos_pthread_mutex_init pthread_mutex_init
#define mos_sem_init sem_init
#define mos_pthread_mutex_lock pthread_mutex_lock
#define mos_pthread_mutex_unlock pthread_mutex_unlock
#define mos_sem_wait sem_wait
#define mos_sem_post sem_post
#define mos_pthread_join pthread_join
#define mos_pthread_detach pthread_detach
#define mos_pthread_self pthread_self
#define mos_pthread_create pthread_create
#endif // _WIN32
