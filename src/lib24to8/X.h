#include <stdint.h>

#undef False
#undef True
#define XLIB_ILLEGAL_ACCESS  1
#include <X11/Xos.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xresource.h>
#include <X11/Xproto.h>
#include <X11/Xatom.h>
#include <X11/cursorfont.h>
#include <X11/keysym.h>
#ifdef HasShape
#include <X11/extensions/shape.h>
#endif
#undef index
#ifdef hpux
#undef SYSV
#endif
#if defined(_AIX) || defined(__hpux)
#define XFD_SET  int
#else
#define XFD_SET  fd_set
#endif

/*
  Default colors declarations.
*/
#define BackgroundColor  "#ccc"  /* gray */
#define BorderColor  "#000"  /* black */
#define ForegroundColor  "#000"  /* black */
#define Pen0Color  "#000"  /* black */
#define Pen1Color  "#000"  /* black */
#define Pen2Color  "#00f"  /* blue */
#define Pen3Color  "#0ff"  /* cyan */
#define Pen4Color  "#0f0"  /* green */
#define Pen5Color  "#ccc"  /* gray */
#define Pen6Color  "#f00"  /* red */
#define Pen7Color  "#f0f"  /* magenta */
#define Pen8Color  "#ff0"  /* yellow */
#define Pen9Color  "#ccc"  /* gray */
/*
  Colormap declarations.
*/
#define UndefinedColormap  0
#define PrivateColormap  1
#define SharedColormap  2
/*
  Define declarations.
*/
#define EditorCommand  "xterm -title 'Image Comment Editor' -e vi %s"
#define MaxNumberFonts  10
#define MaxNumberPens  10
#define PrintCommand  "lpr"
#define RGBColorDatabase  "/usr/lib/X11/rgb.txt"
#define SuspendTime  40
#define XStandardPixel(map,color,dx)  (uint32_t) (map->base_pixel+  \
  ((color.red*map->red_max+(1 << (dx-1)))/((1 << dx)-1))*map->red_mult+  \
  ((color.green*map->green_max+(1 << (dx-1)))/((1 << dx)-1))*map->green_mult+  \
  ((color.blue*map->blue_max+(1 << (dx-1)))/((1 << dx)-1))*map->blue_mult)

/*
  Typedef declarations.
*/
typedef struct _DiversityPacket
{
  unsigned char
    red,
    green,
    blue;

  unsigned short
    index;

  uint32_t 
    count;
} DiversityPacket;

typedef struct _XAnnotateInfo
{
  int
    x,
    y;

  unsigned int
    width,
    height;

  double
    degrees;

  XFontStruct
    *font_info;

  char
    *text,
    geometry[MaxTextLength];

  struct _XAnnotateInfo
    *previous,
    *next;
} XAnnotateInfo;

typedef struct _XPixelInfo
{
  unsigned int
    colors;

  uint32_t 
    *pixels;

  XColor
    foreground_color,
    background_color,
    border_color,
    matte_color,
    highlight_color,
    shadow_color,
    depth_color,
    trough_color,
    pen_color[MaxNumberPens],
    annotate_color;

  unsigned short
    background_index,
    annotate_index;

  GC
   annotate_context,
   highlight_context,
   widget_context;
} XPixelInfo;

typedef struct _XResourceInfo
{
  XrmDatabase
    resource_database;

  ImageInfo
    *image_info;

  unsigned int
    backdrop;

  char
    *background_color,
    *border_color;

  unsigned int
    border_width,
    colormap,
    colorspace,
    debug,
    delay,
    dither;

  char
    *editor_command,
    *font,
    *font_name[MaxNumberFonts],
    *foreground_color;

  int
    gravity;

  char
    *icon_geometry;

  unsigned int
    iconic;

  char
    *image_geometry;

  unsigned int
    magnify;

  char
    *map_type,
    *matte_color;

  unsigned int
    monochrome;

  char
    *name;

  unsigned int
    number_colors;

  char
    *pen_color[MaxNumberPens],
    *print_command;

  char
    *server_name,
    *title;

  unsigned int
    tree_depth,
    update,
    use_pixmap;

  char
    *visual_type,
    *window_id,
    *write_filename;
} XResourceInfo;

typedef struct _XWindowInfo
{
  Window
    id;

  int
    screen;

  Visual
    *visual;

  int
    class,
    depth;

  XVisualInfo
    *visual_info;

  XStandardColormap
    *map_info;

  XPixelInfo
    *pixel_info;

  XFontStruct
    *font_info;

  GC
    annotate_context,
    highlight_context,
    widget_context;

  Cursor
    cursor,
    busy_cursor;

  char
    *name,
    *geometry,
    *icon_name,
    *icon_geometry,
    *clip_geometry;

  uint32_t 
    flags;

  int
    x,
    y;

  unsigned int
    width,
    height,
    min_width,
    min_height,
    width_inc,
    height_inc,
    border_width,
    immutable,
    data;

  XImage
    *ximage,
    *matte_image;

  Pixmap
    highlight_stipple,
    shadow_stipple,
    pixmap,
    matte_pixmap,
    *pixmaps;

  int
    mask;

  XSetWindowAttributes
    attributes;

  XWindowChanges
    window_changes;

  unsigned int
    mapped,
    stasis;
} XWindowInfo;

typedef struct _XWindows
{
  XWindowInfo
    context,
    backdrop,
    icon,
    image,
    info,
    magnify,
    pan,
    command,
    popup;

  Atom
    wm_protocols,
    wm_delete_window,
    wm_take_focus,
    im_protocols,
    im_window_colormap,
    im_former_image,
    im_next_image,
    im_exit;
} XWindows;

/*
  X utilities routines.
*/
extern char
  *XGetResourceClass _Declare((XrmDatabase,char *,char *,char *)),
  *XGetResourceInstance _Declare((XrmDatabase,char *,char *,char *)),
  **XListColors _Declare((char *,int *)),
  *XVisualClassName _Declare((int));

extern Cursor
  XMakeCursor _Declare((Display *,Window,Colormap,char *,char *));

extern Image
  *XGetWindowImage _Declare((Display *,Window,unsigned int,unsigned int)),
  *ReadXImage _Declare((ImageInfo *,unsigned int,unsigned int,unsigned int,
    unsigned int));
       
extern int
  Latin1Compare _Declare((char *,char *)),
  XError _Declare((Display *,XErrorEvent *));

extern unsigned int
  IsTrue _Declare((char *)),
  XAnnotateImage _Declare((Display *,XPixelInfo *,XAnnotateInfo *,unsigned int,
    Image *)),
  XGetWindowColor _Declare((Display *,XColor *)),
  XMakeImage _Declare((Display *,XResourceInfo *,XWindowInfo *,Image *,
    unsigned int,unsigned int)),
  XMakePixmap _Declare((Display *,XResourceInfo *,XWindowInfo *));

extern void
  XBestIconSize _Declare((Display *,XWindowInfo *,Image *)),
  XBestPixel _Declare((Display *,Colormap,XColor *,unsigned int,XColor *)),
  XCheckRefreshWindow _Declare((Display *,XWindowInfo *)),
  XClientMessage _Declare((Display *,Window,Atom,Atom,Time)),
  XDelay _Declare((Display *,uint32_t)),
  XDestroyWindowColors _Declare((Display *,Window)),
  XDisplayInfoString _Declare((Display *,XWindowInfo *,char *)),
  XFreeResources _Declare((Display *,XVisualInfo *,XStandardColormap *,
    XPixelInfo *,XFontStruct *,XResourceInfo *,XWindowInfo *)),
  XFreeStandardColormap _Declare((Display *,XVisualInfo *,XStandardColormap *,
    XPixelInfo *)),
  XGetAnnotateInfo _Declare((XAnnotateInfo *)),
  XGetMapInfo _Declare((XVisualInfo *,Colormap,XStandardColormap *)),
  XGetPixelInfo _Declare((Display *,XVisualInfo *,XStandardColormap *,
    XResourceInfo *,Image *,XPixelInfo *)),
  XGetResourceInfo _Declare((XrmDatabase,char *,XResourceInfo *)),
  XGetWindowInfo _Declare((Display *,XVisualInfo *,XStandardColormap *,
    XPixelInfo *,XFontStruct *,XResourceInfo *,XWindowInfo *)),
  XHighlightLine _Declare((Display *,Window,GC,XSegment *)),
  XHighlightRegion _Declare((Display *,Window,GC,RectangleInfo *)),
  XMakeMagnifyImage _Declare((Display *,XWindows *)),
  XMakeStandardColormap _Declare((Display *,XVisualInfo *,XResourceInfo *,
    Image *,XStandardColormap *,XPixelInfo *)),
  XMakeWindow _Declare((Display *,Window,char **,int,XClassHint *,XWMHints *,
    XWindowInfo *)),
  XRetainWindowColors _Declare((Display *,Window)),
  XSetWindowExtents _Declare((Display *,XWindowInfo *,char *)),
  XRefreshWindow _Declare((Display *,XWindowInfo *,XEvent *));

extern Window
  XClientWindow _Declare((Display *,Window)),
  XSelectWindow _Declare((Display *,RectangleInfo *)),
  XWindowByID _Declare((Display *,Window,uint32_t)),
  XWindowByName _Declare((Display *,Window,char *));

extern XFontStruct
  *XBestFont _Declare((Display *,XResourceInfo *));

extern XVisualInfo
  *XBestVisualInfo _Declare((Display *,XStandardColormap *,XResourceInfo *));

/*
  Variable declarations
*/
extern char
  *client_name;

/*
  Invoke pre-X11R6 ICCCM routines if XlibSpecificationRelease is not 6.
*/
#if XlibSpecificationRelease < 6
#define PRE_R6_ICCCM
#endif
/*
  Invoke pre-X11R5 ICCCM routines if XlibSpecificationRelease is not defined.
*/
#ifndef XlibSpecificationRelease
#define PRE_R5_ICCCM
#endif
/*
  Invoke pre-X11R4 ICCCM routines if PWinGravity is not defined.
*/
#ifndef PWinGravity
#define PRE_R4_ICCCM
#endif
#include "PreRvIcccm.h"
