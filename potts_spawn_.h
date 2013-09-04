/*******************************************************************************/
/*** CREATION FUNCTIONS ***/

  void  readCells();
  void  readCollagen();
  void  putCells();
  void  putCellsHelper(int, int, int);
  void  calculatePerimeter(int);
  void  putCollagen();
  void  putCollagenHelper(int, int, double);

/*******************************************************************************/
/*** Reads us in some cells ***/

void readCells()
{
  int x,y,t;
  FILE* inFile=fopen("lattice_.dat","r");
  if(inFile==NULL){
    printf("\n Nope... missing lattice_.dat\n\n");
    exit(0);
  }
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      fscanf(inFile,"%d",&lattice[i][j][0]);
      if(lattice[i][j][0]>numCells)
        lattice[i][j][0]-=COLLAGEN_OFFSET+PERIMETER_OFFSET*numCells;
      if(lattice[i][j][0]>0){
        cellVolumeList[lattice[i][j][0]][cellVolume[lattice[i][j][0]]][0]=i;
        cellVolumeList[lattice[i][j][0]][cellVolume[lattice[i][j][0]]][1]=j;
        cellVolume[lattice[i][j][0]]++;
      }
    }
  }
  fclose(inFile);
  for(int cell=1;cell<=numCells; cell++)
    calculatePerimeter(cell);
}

/*******************************************************************************/
/*** Reads us in some collagen ***/

void readCollagen()
{
  int x,y,t;
  FILE* inFile=fopen("collagen_.dat","r");
  if(inFile==NULL){
    printf("\n Nope... missing collagen_.dat\n\n");
    exit(0);
  }
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      fscanf(inFile,"%d",&lattice[i][j][1]);
  fclose(inFile);
}

/*******************************************************************************/
/*** Populates an empty upper lattice ***/

void putCells()
{
  if(numCells==1)
    putCellsHelper(N/2,N/2,1);
  else
    for(int cell=1; cell<=numCells; cell++)
      putCellsHelper(rand()%N,rand()%N,cell);
  for(int cell=1;cell<=numCells; cell++)
    calculatePerimeter(cell);
}

/*******************************************************************************/
/*** Places a single cell into the upper lattice ***/

void putCellsHelper(int x0, int y0, int cell)
{
  cellVolume[cell]=0;
  for(int x=-(int)(cellSpawn+0.5); x<=(int)(cellSpawn+0.5); x++){
    int ymax=(int)(sqrt(cellSpawn*cellSpawn-(double)(x*x)+0.5));
    for(int y=-ymax; y<=ymax; y++){
      if(lattice[(N+x0+x)%N][(N+y0+y)%N][0]==0){
        lattice[(N+x0+x)%N][(N+y0+y)%N][0]=cell;
        cellVolumeList[cell][cellVolume[cell]][0]=(N+x0+x)%N;
        cellVolumeList[cell][cellVolume[cell]][1]=(N+y0+y)%N;
        cellVolume[cell]++;
      }
    }
  }
}

/*******************************************************************************/
/*** Calculates cell perimeters ***/

void calculatePerimeter(int cell)
{
  cellPerimeter[cell]=0;
  for(int v=0;v<cellVolume[cell];v++){
    int i=cellVolumeList[cell][v][0];
    int j=cellVolumeList[cell][v][1];
    if( lattice[i][j][0]!=lattice[(i+1)%N][j][0] ||
        lattice[i][j][0]!=lattice[i][(j+1)%N][0] || 
        lattice[i][j][0]!=lattice[(N+i-1)%N][j][0] || 
        lattice[i][j][0]!=lattice[i][(N+j-1)%N][0] )
    {
      cellPerimeterList[cell][cellPerimeter[cell]][0]=i;
      cellPerimeterList[cell][cellPerimeter[cell]][1]=j;
      cellPerimeter[cell]++;
    }
  }
}

/*******************************************************************************/
/*** Puts collagen into the lower lattice ***/

void putCollagen()
{
  for(int n=0; n<numCollagen; n++){
    int x = rand()%N;
    int y = rand()%N;
    double m = tan((double)rand()/((double)RAND_MAX*2.0/3.14159));
    putCollagenHelper(x,y,m);
  }
}

/*******************************************************************************/
/*** Puts a single line of collagen into the lower lattice ***/

void putCollagenHelper(int x0, int y0, double slope)
{
  int sx = 2*(rand()%2)-1;
  int sy = 2*(rand()%2)-1;
  for(int i=0;i<collagenWidth;i++){
    int x  = x0+i*sy;
    int y  = y0-i*sx;
    double error = slope;
    int j=0;
    while(j<N){
      lattice[x][y][1]=1;
      while(error>0.5){
        y=(N+y+sy)%N;
        lattice[x][y][1]=1;
        error-=1.0;
        if(slope>1.0)
          j++;
      }
      x=(N+x+sx)%N;
      error+=slope;
      if(slope<=1.0)
        j++;
    }
  }
}

/*******************************************************************************/
































