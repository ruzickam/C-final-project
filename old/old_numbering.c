// numbering - various shift for various angles
    if ( numbering_angle*57.29578 == 0 )
    {
      x1 = circle_x + sin(numbering_angle)*(circle_r+25);
      y1 = circle_y + cos(numbering_angle)*(circle_r+25);
    }
    if ( numbering_angle*57.29578 < 45 && numbering_angle*57.29578 != 0 )	
    {
      x1 = circle_x + sin(numbering_angle-angle*0.25*(PI/180.0))*(circle_r+23);
      y1 = circle_y + cos(numbering_angle-angle*0.25*(PI/180.0))*(circle_r+23);
    }
    if ( numbering_angle*57.29578 >= 45 && numbering_angle*57.29578 < 90 )
    {
      x1 = circle_x + sin(numbering_angle)*(circle_r+25);
      y1 = circle_y + cos(numbering_angle)*(circle_r+25);
    }
    if ( numbering_angle*57.29578 >= 90 && numbering_angle*57.29578 < 135 )
    {
      x1 = circle_x + sin(numbering_angle+angle*0.27*(PI/180.0))*(circle_r+27);
      y1 = circle_y + cos(numbering_angle+angle*0.27*(PI/180.0))*(circle_r+27);
    }
    if ( numbering_angle*57.29578 >= 135 && numbering_angle*57.29578 < 180 )
    {
      x1 = circle_x + sin(numbering_angle+angle*0.35*(PI/180.0))*(circle_r+32);
      y1 = circle_y + cos(numbering_angle+angle*0.35*(PI/180.0))*(circle_r+32);
    }
    if ( numbering_angle*57.29578 >= 180 && numbering_angle*57.29578 < 225 )
    {
      x1 = circle_x + sin(numbering_angle+angle*0.30*(PI/180.0))*(circle_r+42);
      y1 = circle_y + cos(numbering_angle+angle*0.30*(PI/180.0))*(circle_r+42);
    }
    if ( numbering_angle*57.29578 >= 225 && numbering_angle*57.29578 < 270 )
    {
      x1 = circle_x + sin(numbering_angle+angle*0.27*(PI/180.0))*(circle_r+40);
      y1 = circle_y + cos(numbering_angle+angle*0.27*(PI/180.0))*(circle_r+40);
    }
    if ( numbering_angle*57.29578 >= 270 && numbering_angle*57.29578 < 315 )
    {
      x1 = circle_x + sin(numbering_angle+angle*0.05*(PI/180.0))*(circle_r+35);
      y1 = circle_y + cos(numbering_angle+angle*0.05*(PI/180.0))*(circle_r+35);
    }
    if ( numbering_angle*57.29578 >= 315 && numbering_angle*57.29578 < 360 )
    {
      x1 = circle_x + sin(numbering_angle-angle*0.20*(PI/180.0))*(circle_r+27);
      y1 = circle_y + cos(numbering_angle-angle*0.20*(PI/180.0))*(circle_r+27);
    }