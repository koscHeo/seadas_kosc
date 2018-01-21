#ifndef MAPPROTO_H
#define MAPPROTO_H

int32 get_l3m(char *l3m_path, int16 *syear, int16 *sday, int32 *smsec,
                   int16 *eyear, int16 *eday, int32 *emsec, 
		   char *prod_type, char *l3m_name, uint8 *l3m_data, 
		   unsigned char *palette, meta_struct *meta_l3m);

int32 read_meta(int32 sdfid, meta_struct *meta_l3m);

int32 getattrsz(int32 id, char *attr_name, int32 *nt, int32 *count);

int32 rdattr(int32 sdfid, char *attr_name, void *buf);

#endif /* MAPPROTO_H */
