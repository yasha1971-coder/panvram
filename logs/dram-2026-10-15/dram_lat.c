#define _GNU_SOURCE
#include <errno.h>
#include <inttypes.h>
#include <sched.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>

typedef struct Node { struct Node *next; unsigned char pad[56]; } Node;
_Static_assert(sizeof(Node) == 64, "requires 64-bit pointers and 64-byte nodes");
static uint64_t rng;
static Node * volatile sink;
static void die(const char *s) { perror(s); exit(1); }
static uint64_t number(const char *s) {
    char *end; errno = 0;
    unsigned long long n = strtoull(s, &end, 10);
    if (errno || !*s || *end || *s == '-') { fprintf(stderr,"invalid number: %s\n",s); exit(2); }
    return n;
}
static uint64_t random64(void) {
    uint64_t z = (rng += UINT64_C(0x9e3779b97f4a7c15));
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}
static uint64_t bounded(uint64_t n) {
    uint64_t x, threshold = (uint64_t)(-n) % n;
    do { x = random64(); } while (x < threshold);
    return x % n;
}
static void cycle(Node *a, size_t n) {
    for (size_t i=0; i<n; ++i) a[i].next = a+i;
    /* Sattolo: uniform random single-cycle permutation, not disjoint cycles. */
    for (size_t i=n-1; i>0; --i) {
        size_t j = (size_t)bounded(i);
        Node *tmp=a[i].next; a[i].next=a[j].next; a[j].next=tmp;
    }
}
static uint64_t ns(void) {
    struct timespec t;
    if (clock_gettime(CLOCK_MONOTONIC_RAW,&t)) die("clock_gettime");
    return (uint64_t)t.tv_sec*UINT64_C(1000000000)+(uint64_t)t.tv_nsec;
}
__attribute__((noinline)) static Node *chase(Node *p, uint64_t count) {
    /* Every load depends on the preceding result; volatile prevents dead-code removal. */
    while (count >= 8) {
        p = *(Node * volatile *)p; p = *(Node * volatile *)p;
        p = *(Node * volatile *)p; p = *(Node * volatile *)p;
        p = *(Node * volatile *)p; p = *(Node * volatile *)p;
        p = *(Node * volatile *)p; p = *(Node * volatile *)p;
        count -= 8;
    }
    while (count--) p = *(Node * volatile *)p;
    sink=p; return p;
}
static int compare(const void *a,const void *b) {
    uint64_t x=*(const uint64_t *)a,y=*(const uint64_t *)b;
    return (x>y)-(x<y);
}
static void self_test(void) {
    for (size_t n=2;n<=1024;n*=2) {
        Node *a=calloc(n,sizeof(*a)); unsigned char *seen=calloc(n,1);
        if (!a || !seen) die("calloc");
        for (uint64_t seed=0;seed<20;++seed) {
            rng=seed;cycle(a,n);memset(seen,0,n);Node *p=a;
            for (size_t i=0;i<n;++i) {
                if (p<a || p>=a+n || seen[p-a]) { fputs("cycle test FAIL\n",stderr);exit(1); }
                seen[p-a]=1;p=p->next;
            }
            if (p!=a || chase(a,n)!=a || chase(a,n+1)!=a->next) {
                fputs("chase test FAIL\n",stderr);exit(1);
            }
        }
        free(seen);free(a);
    }
    puts("SELF_TEST PASS cycles=200 (10 sizes x 20 seeds)");
}
int main(int argc,char **argv) {
    int cpu=-1;uint64_t max=UINT64_C(4294967296),minimum=UINT64_C(4194304),seed=20261012;
    for (int i=1;i<argc;++i) {
        if (!strcmp(argv[i],"--self-test")) {self_test();return 0;}
        if (i+1>=argc) { fputs("missing argument\n",stderr);return 2; }
        uint64_t v=number(argv[++i]);
        if (!strcmp(argv[i-1],"--cpu")) {if(v>=CPU_SETSIZE)return 2;cpu=(int)v;}
        else if (!strcmp(argv[i-1],"--max-bytes")) max=v;
        else if (!strcmp(argv[i-1],"--min-accesses")) minimum=v;
        else if (!strcmp(argv[i-1],"--seed")) seed=v;
        else { fputs("unknown option\n",stderr);return 2; }
    }
    if (max<32768 || max>UINT64_C(4294967296) || (max&(max-1)) || !minimum || minimum>UINT64_C(1000000000)) return 2;
    cpu_set_t allowed,one;
    if (sched_getaffinity(0,sizeof(allowed),&allowed))die("sched_getaffinity");
    if (cpu<0) for(int i=0;i<CPU_SETSIZE;++i)if(CPU_ISSET(i,&allowed)){cpu=i;break;}
    if(cpu<0 || !CPU_ISSET(cpu,&allowed)){fputs("CPU outside allowed affinity\n",stderr);return 2;}
    CPU_ZERO(&one);CPU_SET(cpu,&one);
    if(sched_setaffinity(0,sizeof(one),&one))die("sched_setaffinity");
    printf("bytes\tnodes\tcpu\tseed\taccesses_per_trial\twarmups\tmeasurements\tmedian_ns_per_access");
    for(int i=1;i<=9;++i)printf("\ttrial%d_ns_per_access",i);
    puts("\tminor_faults_timed\tmajor_faults_timed");fflush(stdout);
    for(uint64_t bytes=32768;bytes<=max;bytes*=2) {
        size_t n=(size_t)(bytes/64);
        Node *a=mmap(NULL,(size_t)bytes,PROT_READ|PROT_WRITE,MAP_PRIVATE|MAP_ANONYMOUS,-1,0);
        if(a==MAP_FAILED)die("mmap");
        if(madvise(a,(size_t)bytes,MADV_NOHUGEPAGE))die("MADV_NOHUGEPAGE");
        rng=seed^bytes;cycle(a,n); /* first touch after affinity; no auxiliary 4-GiB permutation */
        uint64_t count=minimum>n?minimum:n;
        Node *p=a;for(int i=0;i<3;++i)p=chase(p,count);
        struct rusage before,after;
        if(getrusage(RUSAGE_SELF,&before))die("getrusage");
        uint64_t times[9],sorted[9];
        for(int i=0;i<9;++i) {
            if(sched_getcpu()!=cpu){fputs("affinity changed\n",stderr);return 1;}
            uint64_t start=ns();p=chase(p,count);times[i]=ns()-start;
            if(!times[i]){fputs("zero interval\n",stderr);return 1;}
        }
        if(getrusage(RUSAGE_SELF,&after))die("getrusage");
        if(after.ru_majflt!=before.ru_majflt){fputs("major fault during timing: refused\n",stderr);return 1;}
        memcpy(sorted,times,sizeof(times));qsort(sorted,9,sizeof(*sorted),compare);
        printf("%"PRIu64"\t%zu\t%d\t%"PRIu64"\t%"PRIu64"\t3\t9\t%.6f",bytes,n,cpu,seed,count,(double)sorted[4]/(double)count);
        for(int i=0;i<9;++i)printf("\t%.6f",(double)times[i]/(double)count);
        printf("\t%ld\t%ld\n",after.ru_minflt-before.ru_minflt,after.ru_majflt-before.ru_majflt);fflush(stdout);
        fprintf(stderr,"completed bytes=%"PRIu64" cpu=%d\n",bytes,cpu);
        if(munmap(a,(size_t)bytes))die("munmap");
    }
    return 0;
}
