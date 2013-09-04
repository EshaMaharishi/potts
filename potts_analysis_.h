/*******************************************************************************/
/*** CELL ANALYSIS FUNCTIONS ***/

  void    measureCells();
  double  measureAnisotropy(int);

/*******************************************************************************/
/*** Measure cell properties ***/

void measureCells()
{

  double cellAnisotropy[numCells+1];
  for(int cell=1;cell<=numCells;cell++)
    cellAnisotropy[cell]=measureAnisotropy(cell);

  for(int n=0;n<3;n++){
    avg[n]=0.0;
    dev[n]=0.0;
  }
  for(int cell=1;cell<=numCells;cell++){
    avg[0]+=cellVolume[cell];
    avg[1]+=cellPerimeter[cell];
    avg[2]+=cellAnisotropy[cell];
  }
  for(int n=0;n<3;n++)
    avg[n]/=numCells;
  for(int cell=1;cell<=numCells;cell++){
    dev[0]+=(cellVolume[cell]-avg[0])*(cellVolume[cell]-avg[0]);
    dev[1]+=(cellPerimeter[cell]-avg[1])*(cellPerimeter[cell]-avg[1]);
    dev[2]+=(cellAnisotropy[cell]-avg[2])*(cellAnisotropy[cell]-avg[2]);
  }
  for(int n=0;n<3;n++)
    dev[n]=sqrt(dev[n])/numCells;

  return;
}

/*******************************************************************************/
/*** Measure anisotropy ***/

// NB: THIS WILL FUCK UP IF A CELL GETS BIGGER THAN HALF THE BOX LENGTH.

double measureAnisotropy(int cell)
{

  int    i,j,xi,xf,yi,yf,dx,dy,s;
  int    dist,distmax,distind[4];
  double slope,error;

  /*
    Go through every pair of perimeter spins, find the pair
    which are furthest apart, and record their positions.
  */

  distmax=0;
  for(i=0; i<cellPerimeter[cell]; i++){
    xi=cellPerimeterList[cell][i][0];
    yi=cellPerimeterList[cell][i][1];
    for(j=0; j<cellPerimeter[cell]; j++){
      xf=cellPerimeterList[cell][j][0];
      yf=cellPerimeterList[cell][j][1];
      dx=xf-xi-N*(int)floor((float)(xf-xi)/(float)N+0.499);
      dy=yf-yi-N*(int)floor((float)(yf-yi)/(float)N+0.499);
      dist=dx*dx+dy*dy;
      if(dist>distmax){
        distmax=dist;
        distind[0]=xi;
        distind[1]=xf;
        distind[2]=yi;
        distind[3]=yf;
      }
    }
  }

  xi=distind[0];
  xf=distind[1];
  yi=distind[2];
  yf=distind[3];

  /*
     Find the midpoint.
  */

  dx = (xf-xi)-N*(int)floor((float)(xf-xi)/(float)N+0.499);
  dy = (yf-yi)-N*(int)floor((float)(yf-yi)/(float)N+0.499);

  int xm=(N+xi+dx/2)%N;
  int ym=(N+yi+dy/2)%N;

  /*
    Calculate the slope of the line perpendicular to
    the line between our two perimeter spins.
  */

  if(dy!=0)
    slope = -1.0*(double)dx/(double)dy;
  else
    slope = 2.0*(double)N;

  if(slope>=0.0)
    s=1;
  else
    s=-1;

  /*
    Start at the midpoint, and go out along the perpendicular
    until you find something not in the cell.
  */

  xi = xm;
  yi = ym;
  error = fabs(slope);

  if(lattice[xi][yi][0]!=cell)
    goto donei;

  while(1){
    if(lattice[xi][yi][0]!=cell){
      xi=(N+xi-1)%N;
      goto donei;
    }
    while(error>0.5){
      yi=(N+yi+s)%N;
      error-=1.0;
      if(lattice[xi][yi][0]!=cell){
        yi=(N+yi-s)%N;
        goto donei;
      }
    }
    xi=(xi+1)%N;
    error+=fabs(slope);
  }

  donei:

  xf = xm;
  yf = ym;
  error = fabs(slope);

  if(lattice[xf][yf][0]!=cell)
    goto donef;

  while(1){
    if(lattice[xf][yf][0]!=cell){
      xf=(xf+1)%N;
      goto donef;
    }
    while(error>0.5){
      yf=(N+yf-s)%N;
      error-=1.0;
      if(lattice[xf][yf][0]!=cell){
        yf=(N+yf+s)%N;
        goto donef;
      }
    }
    xf=(N+xf-1)%N;
    error+=fabs(slope);
  }

  donef:

  // Find the length of the perpendicular line.

  dx = xf-xi-N*(int)floor((float)(xf-xi)/(float)N+0.5);
  dy = yf-yi-N*(int)floor((float)(yf-yi)/(float)N+0.5);

  // Return

  double ani;

  if(dx==0 && dy==0)
    ani = (double)distmax;
  else
    ani = sqrt((double)distmax/(double)(dx*dx+dy*dy));

  return ani;
}










































