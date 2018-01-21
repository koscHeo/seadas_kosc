#ifdef PRE_R6_ICCCM
/*
  Compatability defines for pre X11R5 ICCCM.
*/
extern Status
  XInitImage();
#endif

#ifdef PRE_R5_ICCCM
/*
  Compatability defines for pre X11R5 ICCCM.
*/
extern XrmDatabase
  XrmGetDatabase();
#endif

#ifdef PRE_R4_ICCCM
/*
  Compatability defines for pre X11R4 ICCCM.
*/
#ifdef vms
#define XMaxRequestSize(display)  16384
#endif

#define WithdrawnState  0
#define XInductColormap(display,colormap)  XInstallColormap(display,colormap)
#define XUninductColormap(display,colormap)  XUninstallColormap(display,colormap)

typedef struct _XTextProperty
{
  unsigned char
    *value;

  Atom
    encoding;

  int
    format;

  uint32_t
    nitems;
} XTextProperty;

/*
  Pre R4 ICCCM compatibility routines.
*/
char 
  *XResourceManagerString();

extern int
  XWMGeometry();

extern Status
  XGetRGBColormaps(),
  XGetWMName(),
  XReconfigureWMWindow(),
  XSetWMProtocols(),
  XWithdrawWindow();

extern XClassHint
  *XAllocClassHint();

extern XIconSize
  *XAllocIconSize();

extern XSizeHints
  *XAllocSizeHints();

extern XStandardColormap
  *XAllocStandardColormap();

extern XWMHints
  *XAllocWMHints();

extern VisualID
  XVisualIDFromVisual();

extern void
  XrmDestroyDatabase(),
  XSetWMProperties();
#else
#define XInductColormap(display,colormap)
#define XUninductColormap(display,colormap)
#endif
