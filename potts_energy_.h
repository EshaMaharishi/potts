/*******************************************************************************/
/*** HAMILTONIAN FUNCTIONS ***/

  double  inplaneEnergy(int,int);
  double  outplaneEnergy(int,int);
  double  Hamiltonian();
  double  interactionEnergy(int);
  double  volumeEnergy(int);
  double  perimeterEnergy(int);
  double  blobularEnergy(int);

/*******************************************************************************/
/*** Returns the in-plane interaction energy of a lattice site ****/

double inplaneEnergy(int a, int b)
{
  double energy = 0.0;

  if( lattice[a][b][0] != 0 ){
    if( lattice[a][b][0] != lattice[(a+1)%N][b][0] ){
      if( lattice[(a+1)%N][b][0] > 0 )
        return J_cel;
      else
        energy = J_air;
    }
    if( lattice[a][b][0] != lattice[(N+a-1)%N][b][0] ){
      if( lattice[(N+a-1)%N][b][0] > 0 )
        return J_cel;
      else
        energy = J_air;
    }
    if( lattice[a][b][0] != lattice[a][(b+1)%N][0] ){
      if( lattice[a][(b+1)%N][0] > 0 )
        return J_cel;
      else
        energy = J_air;
    }
    if( lattice[a][b][0] != lattice[a][(N+b-1)%N][0] ){
      if( lattice[a][(N+b-1)%N][0] > 0 )
        return J_cel;
      else
        return J_air;
    }
  }

  return energy;
}

/*******************************************************************************/
/*** Returns the out-of-plane interaction energy of a lattice site ****/

double outplaneEnergy(int a, int b)
{
  if( lattice[a][b][0]!=0 && lattice[a][b][1]!=0 )
    return J_col;
  else
    return 0.0;
}

/*******************************************************************************/
/*** Returns the full hamiltonian of the lattice ***/

double Hamiltonian()
{
  double energy = 0.0;
  for(int cell=1;cell<=numCells;cell++){
    energy += interactionEnergy(cell);
    energy += volumeEnergy(cell);
    //energy += perimeterEnergy(cell);
    energy += blobularEnergy(cell);
  }
  return energy;
}

/*******************************************************************************/
/*** Returns the interaction energy of a cell ****/

double interactionEnergy(int cell)
{
  double energy = 0.0;

  for(int p=0;p<cellPerimeter[cell];p++)
    energy+=inplaneEnergy(cellPerimeterList[cell][p][0],cellPerimeterList[cell][p][1]);

  for(int v=0;v<cellVolume[cell];v++)
    energy+=outplaneEnergy(cellVolumeList[cell][v][0],cellVolumeList[cell][v][1]);

  return energy;
}

/*******************************************************************************/
/*** Returns the volume energy of a cell ***/

double volumeEnergy(int cell)
{
  return L_vol*((double)cellVolume[cell]-targetVolume)*((double)cellVolume[cell]-targetVolume);
}

/*******************************************************************************/
/*** Returns the perimeter energy of a cell ***/

double perimeterEnergy(int cell)
{
  return L_per*(double)cellPerimeter[cell]/(double)cellVolume[cell];
}

/*******************************************************************************/
/*** Returns the blobular energy of a cell ***/

double blobularEnergy(int cell)
{

  int energy = 0;
  int number = 0;

  for(int a=0; a<cellPerimeter[cell]; a++){
    for(int b=0; b<cellPerimeter[cell];b++){

      int ai=cellPerimeterList[cell][a][0];
      int aj=cellPerimeterList[cell][a][1];
      int bi=cellPerimeterList[cell][b][0];
      int bj=cellPerimeterList[cell][b][1];

      int dx = bi-ai-N*(int)floor((float)(bi-ai)/(float)N+0.5);
      int sx = (dx>0)-(dx<0);
      int dy = bj-aj-N*(int)floor((float)(bj-aj)/(float)N+0.5);
      int sy = (dy>0)-(dy<0);

      double slope;
      if(dx!=0)
        slope = (double)dy/(double)dx;
      else
        slope = (double)N;

      int x = ai;
      int y = aj;
      double error = fabs(slope);

      do{
        number++;
        if(lattice[x][y][0]!=cell){
          energy++;
          goto done;
        }
        while(error>0.5){
          y=(y+sy+N)%N;
          number++;
          if(lattice[x][y][0]!=cell){
            energy++;
            goto done;
          }
          error=error-1.0;
        }
        x=(x+sx+N)%N;
        error+=fabs(slope);
      }while(x!=(ai+dx+N)%N);

      done:;

    }
  }

  return L_blb * (double)energy / (double)(cellPerimeter[cell]-1);

}

/*******************************************************************************/


























