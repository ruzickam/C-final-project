
// Function for listing the residues on screen + listing of an atom from peptidic backbone
void print_res_backbone()
{
  int i = 0;
  int j = 0;
  char backbone_n[] = " N  ";		// sample atoms from peptidic backbone
  char backbone_ca[] = " CA ";		// --
  char backbone_c[] = " C  ";		// --
  char backbone_o[] = " O  ";		// --
  int backbone_n_number = 0;		// atom number N, CA, C, O
  int backbone_ca_number = 0;
  int backbone_c_number = 0;
  int backbone_o_number = 0;
  
  printf("\n");
  printf("Residue listing, including atom from peptidic backbone:\n");
  printf("\n");
  
  for(i = 0; i < res_count; i++) {
    
    // find backbone atom
    for(j = residues[i].last_atom; j >= residues[i].first_atom; j--) {
      
      if ( strncmp(atoms[j].atom_name, backbone_n, 4) == 0 ) {
	backbone_n_number = atoms[j].atom_number;
      }
      if ( strncmp(atoms[j].atom_name, backbone_ca, 4) == 0 ) {
	backbone_ca_number = atoms[j].atom_number;
      }
      if ( strncmp(atoms[j].atom_name, backbone_c, 4) == 0 ) {
	backbone_c_number = atoms[j].atom_number;
      }
      if ( strncmp(atoms[j].atom_name, backbone_o, 4) == 0 ) {
	backbone_o_number = atoms[j].atom_number;
      }
    }
    
    printf("Residue: %4d %-3.3s, backbone atoms: %4d N, %4d CA, %4d C, %4d O\n", residues[i].residue_number, residues[i].residue_name, backbone_n_number, backbone_ca_number, backbone_c_number, backbone_o_number);
    
    backbone_n_number = 0;
    backbone_ca_number = 0;
    backbone_c_number = 0;
    backbone_o_number = 0;
  }  
}


// Function for drawing the circle with protein sequence
void draw_circle()
{
  int i = 0;
  int win = 0;				// displayed window
  int segments_number = 0;		// number of circle segments
  double angle = 0;			// angle between the segments
  
  segments_number = res_count;
  
  angle = 360.0 / segments_number;
  
  // open the window
  win = g2_open_X11(600,600);
  
  g2_set_line_width(win, 15);
  
  for(i = 0; i < segments_number; i++) {
    g2_pen(win, g2_ink(win, residue_types[residues[i].residue_type].color_r, residue_types[residues[i].residue_type].color_g, residue_types[residues[i].residue_type].color_b) );
    g2_arc(win, 300, 300, 250, 250, angle*i, angle*(i+1));
  }
  
  getchar();
  
  g2_close(win);
}