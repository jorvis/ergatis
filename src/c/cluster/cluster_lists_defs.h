#define TRUE 1
#define FALSE 0
#define LINELEN 400
#define ACCLEN 50

typedef int bool;

/* List of sequences in one cluster */
typedef struct  HitStruct {
   int               valid;
   char              query[ ACCLEN ];
   char              hit  [ ACCLEN ];
   struct HitStruct  *next;
}  HitStruct;
typedef HitStruct *HITPTR;

HITPTR PAIR_LIST_FIRST, PAIR_LIST_LAST;
