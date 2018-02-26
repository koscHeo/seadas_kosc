#ifndef L1A_SEAWIFS_H
#define L1A_SEAWIFS_H

int openl1a_seawifs(filehandle *file);
int readl1a_seawifs(filehandle *file, int32 recnum, l1str *l1rec);
int readl1a_lonlat_seawifs(filehandle *file, int32 recnum, l1str *l1rec);
int closel1a_seawifs(filehandle *file);


#endif
