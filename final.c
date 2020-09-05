#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <g2.h>
#include <g2_X11.h>
#include <math.h>

/* The final project - course C2160
 * Circular diagram of the distance between the centre of the molecule and residues 
 *  
 * Use the program this way: ./final INPUT_FILE
 * 
 * Author: Michal Ruzicka
 * Published: 2015 
 */

// max. number of atoms
#define MAX_ATOMS 100000

// max. number of residues
#define MAX_RESIDUE 10000

// buffer size for loading from the PDB file
#define BUF_SIZE 1000

// number of residual types for proteins
#define RESD_TYPES_COUNT 21

// PI constant
#define PI 3.14159265358979323846

//struct definitions
typedef struct		// struct ATOM (one ATOM record in the PDB file)
{
  char record_name[7];	// record name for ATOM 
  int atom_number;	// atom number
  char atom_name[5];	// atom name
  char alt_loc;		// alternative atom position
  char residue_name[4];	// residue name
  char chain_id;	// chain ID for proteins
  int residue_number;	// residue number in the sequence
  char i_code;		// code for residue insert
  double x;		// coordinate x in ATOM
  double y;		// coordinate y in ATOM
  double z;		// coordinate z in ATOM
  double occupancy;	// occupancy
  double temp_factor;	// thermal factor
  char element_name[3];	// element name
  char formal_charge[3];// formal charge
} ATOM;

typedef struct		// struct RESIDUE
{
  int first_atom;	// index for the first atom from the residue
  int last_atom;	// index for the last atom from the residue
  char residue_name[4];	// residue name
  int residue_number;	// residue number in the sequence
  int residue_type;	// aminoacid type (0-20)
  int atom_c_alpha;	// C-alpha atom
  double atom_c_alpha_x;// coordinate x for C-alpha atom
  double atom_c_alpha_y;// coordinate y for C-alpha atom
  double atom_c_alpha_z;// coordinate z for C-alpha atom
} RESIDUE;

typedef struct		// struct RESIDUE_TYPE
{
  char code3[4];	// 3-letter abbr
  char code1;		// 1-letter abbr
  double color_r;	// color "red" from RGB
  double color_g;	// color "green" from RGB
  double color_b;	// color "blue" from RGB
} RESIDUE_TYPE;

// definition of global variables
int atom_count = 0;				// counter for ATOM records
int res_count = 0;				// counter for RESIDUE records
char file[150] = "";	 			// path to the input PDB file
ATOM atoms[MAX_ATOMS];				// ATOM records
RESIDUE residues[MAX_RESIDUE];			// RESIDUE records
RESIDUE_TYPE residue_types[RESD_TYPES_COUNT] =	// RESIDUE_TYPE records
{
  {"UNK", 'X', 153/255.0, 153/255.0, 153/255.0},
  {"ALA", 'A', 204/255.0, 255/255.0, 255/255.0},
  {"ARG", 'R', 230/255.0,   6/255.0,   6/255.0},
  {"ASN", 'N', 255/255.0, 153/255.0,   0/255.0},
  {"ASP", 'D', 255/255.0, 204/255.0, 153/255.0},
  {"CYS", 'C',   0/255.0, 255/255.0, 255/255.0},
  {"GLN", 'Q', 255/255.0, 102/255.0,   0/255.0},
  {"GLU", 'E', 255/255.0, 204/255.0,   0/255.0},
  {"GLY", 'G',   0/255.0, 255/255.0,   0/255.0},
  {"HIS", 'H', 255/255.0, 255/255.0, 153/255.0},
  {"ILE", 'I',   0/255.0,   0/255.0, 128/255.0},
  {"LEU", 'L',  51/255.0, 102/255.0, 255/255.0},
  {"LYS", 'K', 198/255.0,   6/255.0,   0/255.0},
  {"MET", 'M', 153/255.0, 204/255.0, 255/255.0},
  {"PHE", 'F',   0/255.0, 204/255.0, 255/255.0},
  {"PRO", 'P', 255/255.0, 255/255.0,   0/255.0},
  {"SER", 'S', 204/255.0, 255/255.0, 153/255.0},
  {"THR", 'T',   0/255.0, 255/255.0, 153/255.0},
  {"TRP", 'W', 204/255.0, 153/255.0, 255/255.0},
  {"TYR", 'Y', 204/255.0, 255/255.0, 204/255.0},
  {"VAL", 'V',   0/255.0,   0/255.0, 255/255.0}
};

// declaration of functions (function main is the first one)
int read_file(char file2read[]);							// Function for reading of the input PDB file, return value 0 = loading was OK, 1 = error during loading
void load_res();									// Function for residue loading
void center_of_gravity(double *p_centerX, double *p_centerY, double *p_centerZ);	// Function for the center of the molecule calculation, the centre of gravity is calculated (all atoms have the same weight)
void get_points_distance(double point1X, double point1Y, double point1Z, double point2X, double point2Y, double point2Z, double *p_distance);	// Function that returns a distance between 2 points
void draw_scheme();									// Function that draws circular diagram with a C-alpha distance from the centre for each atom


int main(int argc, char *argv[])
{
  int a = 0;			// return value for the read_file function
  
  if ( argc < 2 )
  {
    printf("No parameters! Use the program this way: ./final INPUT_FILE\n");
    return 1;
  }
  
  if ( strcmp(argv[1],"-h") == 0 )
  {
    printf("Use the program this way: ./final INPUT_FILE\n");
    return 0;
  }
  
  // loading of the arguments
  strncpy(file, argv[1], strlen(argv[1]));
  file[strlen(argv[1])] = '\0';

  // loading of the PDB file
  a = read_file(file);
  if ( a == 1 )
    return 1;
  
  // loading of the residues
  load_res();
  
  // drawing of the circular diagram with a C-alpha distance from the centre for each atom
  draw_scheme();
  
  return 0;
}


// Function for reading of the input PDB file, return value 0 = loading was OK, 1 = error during loading
int read_file(char file2read[])
{
  FILE *fread = NULL; 		// var FILE that identifies the file for reading
  char buf[BUF_SIZE] = "";	// buffer for line loading
  char s[30] = "";		// string var as a help
  
  fread = fopen(file2read, "r");
  if (fread == NULL)
  {
    printf("Cannot open input file %s!\n", file2read);
    if (errno != 0)
      printf("	Error explanation: %s\n", strerror(errno));
    return 1;
  }
  
  while (feof(fread) == 0)
  {    
    if ( atom_count >= MAX_ATOMS )
    {
      printf("PDB file is too big (not all ATOM records were loaded)!\n");
      break;
    }
    
    memset(buf, '\0', BUF_SIZE);
    if ( fgets(buf, BUF_SIZE, fread ) == NULL )	// loading of one line, in case of a failure the cycle is ended
    {
      break;
    }
    
    if ( strncmp(buf, "ATOM", 4) == 0 )	// loading only of record ATOM
    {
      // record name 
      strncpy(atoms[atom_count].record_name, buf, 6);
      atoms[atom_count].record_name[6] = '\0';      
      
      // atom number
      strncpy(s, buf+6, 5);
      s[5] = '\0';
      sscanf(s, "%d", &atoms[atom_count].atom_number);
      
      // atom name
      strncpy(atoms[atom_count].atom_name, buf+12, 4);
      atoms[atom_count].atom_name[4] = '\0';
      
      //  alternative atom position
      atoms[atom_count].alt_loc = buf[16];
      
      // residue name
      strncpy(atoms[atom_count].residue_name, buf+17, 3);
      atoms[atom_count].residue_name[3] = '\0';
      
      // chain ID for proteins
      atoms[atom_count].chain_id = buf[21];
      
      // residue number in the sequence
      strncpy(s, buf+22, 4);
      s[4] = '\0';
      sscanf(s, "%d", &atoms[atom_count].residue_number);
      
      // code for residue insert
      atoms[atom_count].i_code = buf[26];
      
      // coordinate x in ATOM
      strncpy(s, buf+30, 8);
      s[8] = '\0';
      sscanf(s, "%lf", &atoms[atom_count].x);
      
      // coordinate y in ATOM
      strncpy(s, buf+38, 8);
      s[8] = '\0';
      sscanf(s, "%lf", &atoms[atom_count].y);
      
      // coordinate z in ATOM
      strncpy(s, buf+46, 8);
      s[8] = '\0';
      sscanf(s, "%lf", &atoms[atom_count].z);
      
      // occupancy
      strncpy(s, buf+54, 6);
      s[6] = '\0';
      sscanf(s, "%lf", &atoms[atom_count].occupancy);
      
      // thermal factor
      strncpy(s, buf+60, 6);
      s[6] = '\0';
      sscanf(s, "%lf", &atoms[atom_count].temp_factor);
      
      // element name
      strncpy(atoms[atom_count].element_name, buf+76, 2);
      atoms[atom_count].element_name[2] = '\0';
      
      // formal charge
      strncpy(atoms[atom_count].formal_charge, buf+78, 2);
      atoms[atom_count].formal_charge[2] = '\0';      
      
      atom_count++;
    }    
  }
  
  if ( fclose(fread) == EOF )
  {
    printf("Error closing the file %s\n", file2read);
    if (errno != 0)
      printf("	Error explanation: %s\n", strerror(errno));
    return 1;
  }
  fread = NULL;
  
  return 0;
}  


// Function for residue loading
void load_res()
{
  int i = 0;
  int j = 0;
  int first_atom = 0;		// index for the first atom from the residue
  
  for(i = 0; i < atom_count; i++) {
    if ( atoms[i].residue_number != atoms[i+1].residue_number )  // looking for the interface between two residues
    {
      if ( res_count >= MAX_RESIDUE )
      {
	printf("PDB file is too big (not all RESIDUE records were loaded)!\n");
	break;
      }
      
      residues[res_count].residue_number = atoms[i].residue_number;
      
      strncpy(residues[res_count].residue_name, atoms[i].residue_name, 3);
      residues[res_count].residue_name[3] = '\0';
      
      residues[res_count].last_atom = i;
      
      residues[res_count].first_atom = first_atom;
      first_atom = i+1;
      
      // residue_type loading
      for(j = 0; j < RESD_TYPES_COUNT; j++) {
	if ( strncmp(residues[res_count].residue_name, residue_types[j].code3, 3) == 0 )
	{
	  residues[res_count].residue_type = j;
	}
      }  
      
      // looking for the C-alpha atom in the residue
      for(j = residues[res_count].first_atom; j <= residues[res_count].last_atom; j++) {
	if ( strncmp(atoms[j].atom_name, " CA ", 4) == 0 ) 
	{
	  residues[res_count].atom_c_alpha = j;
	  residues[res_count].atom_c_alpha_x = atoms[j].x;
	  residues[res_count].atom_c_alpha_y = atoms[j].y;
	  residues[res_count].atom_c_alpha_z = atoms[j].z;
	  break;
	}
      }
      
      res_count++;
    }  
  }
}


// Function for the center of the molecule calculation, the centre of gravity is calculated (all atoms have the same weight)
void center_of_gravity(double *p_centerX, double *p_centerY, double *p_centerZ)
{
  int i = 0;
  double x_sum = 0.0;		// sum of x coordinates
  double y_sum = 0.0;		// sum of y coordinates
  double z_sum = 0.0;		// sum of z coordinates
  
  for(i = 0; i < atom_count; i++) {
    x_sum += atoms[i].x;
    y_sum += atoms[i].y;
    z_sum += atoms[i].z;
  }
  
  *p_centerX = x_sum / atom_count;
  *p_centerY = y_sum / atom_count;
  *p_centerZ = z_sum / atom_count;
}


// Function that returns a distance between 2 points
void get_points_distance(double point1X, double point1Y, double point1Z, double point2X, double point2Y, double point2Z, double *p_distance)
{
  *p_distance = sqrt( (point1X - point2X)*(point1X - point2X) + (point1Y - point2Y)*(point1Y - point2Y) + (point1Z - point2Z)*(point1Z - point2Z) );
}


// Function that draws circular diagram with a C-alpha distance from the centre for each atom
void draw_scheme()
{
  int i = 0;
  double centerX = 0.0;				// centre of the molecule (coordinate x)
  double centerY = 0.0;				// centre of the molecule (coordinate y)
  double centerZ = 0.0;				// centre of the molecule (coordinate z)
  double distances[MAX_RESIDUE] = {0.0};	// distance between a residue and the centre
  double max_distance = 0.0;			// max value of distance C-alpha - centre
  int win = 0;					// displayed window
  int win_x = 800;				// x window
  int win_y = 700;				// y window
  int circle_x = 320;				// coordinate x for the circle
  int circle_y = 350;				// coordinate y for the circle
  int circle_r = 250;				// circle radius
  int segments_number = 0;			// number of segments in the circle
  double angle = 0;				// angle between segments
  double x1 = 0, y1 = 0, x2 = 0, y2 = 0;	// coordinates for numbering lines
  int numbering_lines = 0;			// number of numbering lines (numbering every 10 residues)
  double numbering_angle = 0.0;		// angle between numbering lines
  char numbering_label[5] = "";		// numbering label
  
  // calculation of the centre of the molecule coordinates
  center_of_gravity(&centerX, &centerY, &centerZ);
  
  // calculation of the distance from the centre for each residue + calc of max distances
  for(i = 0; i < res_count; i++) {
    get_points_distance(centerX, centerY, centerZ, residues[i].atom_c_alpha_x, residues[i].atom_c_alpha_y, residues[i].atom_c_alpha_z, &distances[i]);
    if ( distances[i] > max_distance )
    {
      max_distance = distances[i];
    } 
    //printf("%s%i: %f\n", residues[i].residue_name, residues[i].residue_number, distances[i]);
  }
  
  // -----------------
  // DRAW SCHEME
  // -----------------
  
  // definition of number of segments and angles between them
  segments_number = res_count;
  angle = 360.0 / segments_number;
  
  // open the window for graphics
  win = g2_open_X11(win_x,win_y);
  
  // draw the circle
  g2_set_line_width(win, 15);
  for(i = 0; i < segments_number; i++) {
    // set up the residue color
    g2_pen(win, g2_ink(win, residue_types[residues[i].residue_type].color_r, residue_types[residues[i].residue_type].color_g, residue_types[residues[i].residue_type].color_b) );
    // draw arc
    g2_arc(win, circle_x, circle_y, circle_r, circle_r, ((90+angle/2)-angle*(i+1)), ((90+angle/2)-angle*i) );
  }
  
  // numbering (every 10th residue)
  g2_pen(win, 1);
  g2_set_line_width(win, 1);
  // calc number of lines
  numbering_lines = (segments_number / 10) + 1;
  for(i = 0; i < numbering_lines; i++) {
    // between the 1st and 10th segment, the angle is smaller
    if ( i == 0 || i == 1 )
    {
      numbering_angle = angle*(PI/180.0)*9.0*i;// factor PI/180 converts degrees to radians
    }
    else
    {
      numbering_angle = angle*(PI/180.0)*(19.0+((i-2)*10.0));
    }
    x1 = circle_x + sin(numbering_angle)*(circle_r+20);
    y1 = circle_y + cos(numbering_angle)*(circle_r+20);
    // line from centre to the edge
    g2_line(win, circle_x, circle_y, x1, y1);
    // numbering
    x1 = circle_x-8 + sin(numbering_angle)*(circle_r+35);
    y1 = circle_y-6 + cos(numbering_angle)*(circle_r+35);
    // calc numbers of the numbering
    if ( i == 0 )
    {
      snprintf(numbering_label, 5, "%i", residues[0].residue_number);
    }
    else
    {
      snprintf(numbering_label, 5, "%i", residues[i*10-1].residue_number);
    }
    g2_string(win, x1, y1, numbering_label);
  }
  
  // line for the distance from residues to the centre
  g2_pen(win, 19);
  g2_set_line_width(win, 1);
  for(i = 0; i < res_count; i++) {
    // coordinate for a residue (x1, y1)
    numbering_angle = angle*(PI/180.0)*i;// faktor PI/180 prevadi stupne na radiany
    x1 = circle_x + sin(numbering_angle)*(circle_r-7.5)*(distances[i]/max_distance);
    y1 = circle_y + cos(numbering_angle)*(circle_r-7.5)*(distances[i]/max_distance);
    
    // draw red junction
    if ( i != 0 ) {
      g2_line(win, x1, y1, x2, y2);
    }
    
    // transfer coordinates for next residue junction
    x2 = x1;
    y2 = y1;
  }
  // draw the last line (red junction between the first and last residue)
  numbering_angle = angle*(PI/180.0)*0;// factor PI/180 converts degrees to radians
  x1 = circle_x + sin(numbering_angle)*(circle_r-7.5)*(distances[0]/max_distance);
  y1 = circle_y + cos(numbering_angle)*(circle_r-7.5)*(distances[0]/max_distance);
  g2_line(win, x1, y1, x2, y2);
  
  // label
  g2_pen(win, 1);
  g2_set_font_size(win, 16);
  g2_string(win, 30, circle_y+circle_r+60, "Distance between a C-alpha from the centre of mass");
  
  // label for PDB file name
  g2_pen(win, 1);
  g2_set_font_size(win, 16);
  g2_string(win, 30, circle_y-circle_r-75, "PDB file:");
  g2_string(win, 120, circle_y-circle_r-75, file);
  
  // key with residues
  g2_set_line_width(win, 0);
  for(i = 0; i < RESD_TYPES_COUNT; i++) {  
    // set up color acoording to the residue type
    g2_pen(win, g2_ink(win, residue_types[i].color_r, residue_types[i].color_g, residue_types[i].color_b) );
    // draw rectangle
    g2_filled_rectangle(win, circle_x+circle_r+100, circle_y+circle_r-(i*25), circle_x+circle_r+125, circle_y+circle_r-(i*25)+15);
    // 3-letter abbr
    g2_pen(win, 1);
    g2_string(win,  circle_x+circle_r+135, circle_y+circle_r-(i*25), residue_types[i].code3);
  }  
  
  getchar();
  
  g2_close(win);
}



