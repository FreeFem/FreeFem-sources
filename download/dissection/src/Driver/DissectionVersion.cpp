void
DissectionVersion (
int * const                 vern,
int * const                 rels,
int * const                 ptch)
{
  *versptr = SCOTCH_VERSION;
  *relaptr = SCOTCH_RELEASE;
  *patcptr = SCOTCH_PATCHLEVEL;
}
