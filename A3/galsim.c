// Copilot was used to help better comment the code
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include <X11/Xlib.h>
#include <X11/keysym.h>
#include <unistd.h>

// Adding the graphics functions

void FlushDisplay(void);
void CloseDisplay(void);
int CheckForQuit(void);
void Refresh(void);
void ClearScreen(void);
void DrawCircle(float x, float y, float W, float H, float radius, float color);
void DrawRectangle(float x, float y, float W, float H, float dx, float dy, float color);
void InitializeGraphics(char *command, int windowWidth, int windowHeight);
void SetCAxes(float cmin, float cmax);


typedef struct {
    double x;
    double y;
} Vector2D;

// Function prototypes
Vector2D get_force_on_body(const int nstars, const int G, const float e0, int i, Vector2D* position, const double* mass);
double get_rij(int i, int j, Vector2D* position);
Vector2D get_position_vector(int i, int j, Vector2D* position);
void update_velocity_and_position(int i, const int stepsize, Vector2D* velocity, Vector2D* position, Vector2D* F, const double* mass);



int main(int argc, char* argv[]) {
    if(argc != 6) {
        fprintf(stderr, "Usage: %s <Number of stars> <input file> <Number of timesteps> <size of timesteps> <graphics>\n", argv[0]);  return 1;
    }
    // Read in command line arguments
    // Make everything const for optimization
    const int nstars = atoi(argv[1]);
    const char* input_file = argv[2];
    const int nsteps = atoi(argv[3]);
    const int stepsize = atoi(argv[4]);
    const bool graphics = atoi(argv[5]); // 0 or 1 for false or true

    const int G = 100/nstars;
    const double e0 = 0.001; // Softening factor 10^-3

    // Create arrays to store the data and Vector2D to store the forces
    Vector2D* position = (Vector2D*) malloc(nstars * sizeof(Vector2D));
    const double* mass = (double*) malloc(nstars * sizeof(double));
    Vector2D* velocity = (Vector2D*) malloc(nstars * sizeof(Vector2D));
    const bool* brightness = (double*) malloc(nstars * sizeof(double));

    if (position == NULL || mass == NULL || velocity == NULL || brightness == NULL) {
        fprintf(stderr, "Error allocating memory\n");
        return 1;
    }

    // Read in the data
    FILE* file = fopen(input_file, "r");
    if(file == NULL) {
        fprintf(stderr, "Error opening file\n");
        return 1;
    }

    // Read in the input file
    // Input file has structure:
    /*particle 0 position x
    particle 0 position y
    particle 0 mass
    particle 0 velocity x
    particle 0 velocity y
    particle 0 brightness
    particle 1 position x*/
    // The input file is binary
    for (int i = 0; i < nstars; i++) {
        fread(&position[i], sizeof(Vector2D), 1, file);
        fread(&mass[i], sizeof(double), 1, file);
        fread(&velocity[i], sizeof(Vector2D), 1, file);
        fread(&brightness[i], sizeof(bool), 1, file);
    }
    fclose(file);

    if (graphics) {
        InitializeGraphics(argv[0], 800, 800);
        SetCAxes(0, 1);
    }

    for (int step = 0; step < nsteps; step++) {
        for (int i = 0; i < nstars; i++) {
            Vector2D F = get_force_on_body(nstars, G, e0, i, position, mass);
            update_velocity_and_position(i, stepsize, velocity, position, &F, mass);
        }
        if (graphics) {
            ClearScreen();
            for (int i = 0; i < nstars; i++) {
                DrawCircle(position[i].x, position[i].y, 1, 1, 0.01, 0.5);
            }
            Refresh();
            usleep(30000); // Sleep for 30ms to control the animation speed
        }
    }

    if (graphics) {
        FlushDisplay();
        CloseDisplay();
    }

    // Output the data in a binary file
    FILE* output = fopen("result.gal", "wb");
    if(output == NULL) {
        fprintf(stderr, "Error opening file\n");
        return 1;
    }
    for (int i = 0; i < nstars; i++) {
        fwrite(&position[i], sizeof(Vector2D), 1, output);
        fwrite(&mass[i], sizeof(double), 1, output);
        fwrite(&velocity[i], sizeof(Vector2D), 1, output);
        fwrite(&brightness[i], sizeof(bool), 1, output);
    }

    // Free the memory
    free(position);
    free(mass);
    free(velocity);
    free(brightness);
    // Jag bytte till test

}
// Function to calculate the force on a body
Vector2D get_force_on_body(const int nstars, const int G, const float e0, int i, Vector2D* position, const double* mass) {
    Vector2D F;
    F.x = -G * mass[i];
    F.y = -G * mass[i]; 
    for (int j = 0; j < nstars; j++) {
        if (i != j) {
            double rij = get_rij(i, j, position);
            Vector2D vector = get_position_vector(i, j, position);
            double temp = mass[j] / pow(rij + e0, 3);
            F.x *= temp * vector.x;
            F.y *= temp * vector.y;
        }
    }
    return F;
}

// Function to calculate the distance between two bodies
double get_rij(int i, int j, Vector2D* position) {
    return sqrt(pow(position[i].x - position[j].x, 2) + pow(position[i].y - position[j].y, 2));
}

// Function to calculate the position vector between two bodies
Vector2D get_position_vector(int i, int j, Vector2D* position) {
    Vector2D rij;
    rij.x = position[i].x - position[j].x;
    rij.y = position[i].y - position[j].y;
    return rij;
}

// Function to calculate the acceleration of a body
void update_velocity_and_position(int i, const int stepsize, Vector2D* velocity, Vector2D* position, Vector2D* F, const double* mass) {
    velocity[i+1].x += stepsize * F[i].x / mass[i];
    velocity[i+1].y += stepsize * F[i].y / mass[i];
    position[i+1].x += stepsize * velocity[i+1].x;
    position[i+1].y += stepsize * velocity[i+1].y;
}

int main_update(int nstars, int G, float e0, int i, Vector2D* position, const double* mass, Vector2D* velocity, int stepsize) {
    Vector2D F = get_force_on_body(nstars, G, e0, i, position, mass);
    update_velocity_and_position(i, stepsize, velocity, position, &F, mass);
    return 0;
}




// Adding all the code for the graphics

#define NUMCOLORS 512

Display *global_display_ptr;

Window win;
Pixmap pixmap;
XEvent report;
GC gc;
unsigned black, white;
unsigned width, height;
unsigned colors[NUMCOLORS];
float caxis[2];

/*
 * function: create_simple_window. Creates a window with a white background
 *           in the given size.
 * input:    display, size of the window (in pixels), and location of the window
 *           (in pixels).
 * output:   the window's ID.
 * notes:    window is created with a black border, 2 pixels wide.
 *           the window is automatically mapped after its creation.
 */
Window create_simple_window(Display* display, int width, int height, int x, int y)
{
  int screen_num = DefaultScreen(display);
  int win_border_width = 2;
  Window win;

  /* create a simple window, as a direct child of the screen's */
  /* root window. Use the screen's black and white colors as   */
  /* the foreground and background colors of the window,       */
  /* respectively. Place the new window's top-left corner at   */
  /* the given 'x,y' coordinates.                              */
  win = XCreateSimpleWindow(display, RootWindow(display, screen_num),
                            x, y, width, height, win_border_width,
                            BlackPixel(display, screen_num),
                            WhitePixel(display, screen_num));

  /* make the window actually appear on the screen. */
  XMapWindow(display, win);

  /* flush all pending requests to the X server. */
  XFlush(display);

  return win;
}

GC create_gc(Display* display, Window win, int reverse_video)
{
  GC gc;				/* handle of newly created GC.  */
  unsigned long valuemask = 0;		/* which values in 'values' to  */
					/* check when creating the GC.  */
  XGCValues values;			/* initial values for the GC.   */
  unsigned int line_width = 1;		/* line width for the GC.       */
  int line_style = LineSolid;		/* style for lines drawing and  */
  int cap_style = CapButt;		/* style of the line's edje and */
  int join_style = JoinBevel;		/*  joined lines.		*/
  int screen_num = DefaultScreen(display);

  gc = XCreateGC(display, win, valuemask, &values);
  if (gc < 0) {
	fprintf(stderr, "XCreateGC: \n");
  }

  /* allocate foreground and background colors for this GC. */
  if (reverse_video) {
    XSetForeground(display, gc, WhitePixel(display, screen_num));
    XSetBackground(display, gc, BlackPixel(display, screen_num));
  }
  else {
    XSetForeground(display, gc, BlackPixel(display, screen_num));
    XSetBackground(display, gc, WhitePixel(display, screen_num));
  }

  /* define the style of lines that will be drawn using this GC. */
  XSetLineAttributes(display, gc,
                     line_width, line_style, cap_style, join_style);

  /* define the fill style for the GC. to be 'solid filling'. */
  XSetFillStyle(display, gc, FillSolid);

  return gc;
}

void InitializeGraphics(char *command, int windowWidth, int windowHeight) {
  int i;
  Screen *screen;
  int screen_num;		/* number of screen to place the window on.  */
  /* height and width of the X display
  unsigned int display_width,
               display_height;
  */
  char *display_name = getenv("DISPLAY");  /* address of the X display.      */
  Colormap screen_colormap;     /* color map to use for allocating colors.   */
  XColor color;

  /* Set the height and width of the window */
  width=windowWidth;
  height=windowHeight;

  /* open connection with the X server. */
  global_display_ptr = XOpenDisplay(display_name);
  if (global_display_ptr == NULL) {
    fprintf(stderr, "%s: cannot connect to X server '%s'\n", command, display_name);
    exit(1);
  }

  /* get the geometry of the default screen for our display. */
  screen = ScreenOfDisplay(global_display_ptr, 0);
  screen_num = XScreenNumberOfScreen(screen);
  /*
  display_width = DisplayWidth(global_display_ptr, screen_num);
  display_height = DisplayHeight(global_display_ptr, screen_num);
  */

  /* create a simple window, as a direct child of the screen's   */
  /* root window. Use the screen's white color as the background */
  /* color of the window. Place the new window's top-left corner */
  /* at the given 'x,y' coordinates.                             */
  win = create_simple_window(global_display_ptr, width, height, 0, 0);
  pixmap = XCreatePixmap(global_display_ptr, win, width, height, DefaultDepthOfScreen(screen));

  /* allocate a new GC (graphics context) for drawing in the window. */
  gc = create_gc(global_display_ptr, win, 0);
  XSync(global_display_ptr, False);

  /* get access to the screen's color map. */
  screen_colormap = DefaultColormap(global_display_ptr, DefaultScreen(global_display_ptr));

  black = BlackPixel(global_display_ptr,screen_num);
  white = WhitePixel(global_display_ptr,screen_num);

  /* Load the jet color map */
  for(i=0;i<NUMCOLORS;i++) {
    color.red=((double)(NUMCOLORS-i)/(double)NUMCOLORS)*0xFFFF;
    color.blue=color.red;
    color.green=color.red;
    XAllocColor(global_display_ptr, screen_colormap,&color);    
    colors[i]=color.pixel;
  }
  SetCAxes(0,1);

  /* Set up the window to wait for a quit signal */
  XSelectInput(global_display_ptr, win, ExposureMask | KeyPressMask | ButtonPressMask | ButtonReleaseMask );
  XCopyArea(global_display_ptr, pixmap, win, gc, 0, 0, width, height, 0, 0);
  XMaskEvent(global_display_ptr, ExposureMask, &report);
}

void SetCAxes(float cmin, float cmax) {
  caxis[0]=cmin;
  caxis[1]=cmax;
}

int CheckForQuit(void) {
  int keysym;

  if(XCheckMaskEvent(global_display_ptr,KeyPressMask,&report)) {
    switch(keysym=XLookupKeysym(&report.xkey, 0)) {
    case XK_q:
      return 1;
      break;
    default:
      break;
    }
  }
  return 0;
}

void Refresh(void) {
  XCopyArea(global_display_ptr, pixmap, win, gc, 0, 0, width, height, 0, 0);
  XFlush(global_display_ptr);
}

void ClearScreen(void) {
  XSetForeground(global_display_ptr, gc, black);
  XFillRectangle(global_display_ptr, pixmap, gc, 0, 0, width, height);
}

void DrawCircle(float x, float y, float W, float H, float radius, float color) {
  int i=(int)((x-radius)/W*width);
  int j=height-(int)((y+radius)/H*height);
  int arcrad=2*(int)(radius/W*width);
  int icolor;

  if(color>=caxis[1])
    icolor=NUMCOLORS-1;
  else if(color<caxis[0])
    icolor=0;
  else
    icolor=(int)((color-caxis[0])/(caxis[1]-caxis[0])*(float)NUMCOLORS);

  XSetForeground(global_display_ptr, gc, colors[icolor]);
  XFillArc(global_display_ptr, pixmap, gc, i, j, arcrad, arcrad, 0, 64*360);
}

void DrawRectangle(float x, float y, float W, float H, float dx, float dy, float color) {
  int i=(int)(x/W*width);
  int j=height-(int)((y+dy)/H*height);
  int w=(int)(dx/W*width);
  int h=(int)(dy/H*height);
  int icolor;

  if(color>=caxis[1])
    icolor=NUMCOLORS-1;
  else if(color<caxis[0])
    icolor=0;
  else
    icolor=(int)((color-caxis[0])/(caxis[1]-caxis[0])*(float)NUMCOLORS);

  XSetForeground(global_display_ptr, gc, colors[icolor]);
  XDrawRectangle(global_display_ptr, pixmap, gc, i, j, w, h);
}

void FlushDisplay(void) {
  XFlush(global_display_ptr);
}

void CloseDisplay(void) {
  XCloseDisplay(global_display_ptr);
}
