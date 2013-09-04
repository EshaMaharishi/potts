/*******************************************************************************/
/*** SPIN FLIP FUNCTIONS ***/

  int    flip();
  void   choose();
  void   adjustVolumes();
  void   adjustPerimeters();
  bool   maintainsContiguity();

  int    iSite;
  int    jSite;
  int    oldCell;
  int    newCell;

/*******************************************************************************/
/*** Flip a spin and accept or reject ***/

int flip()
{

  // Choose a spin

  choose();

  // Old energy

  double deltaEnergy = 0.0;

  deltaEnergy -= outplaneEnergy(iSite,jSite);

  deltaEnergy -= inplaneEnergy(iSite,jSite);
  deltaEnergy -= inplaneEnergy((iSite+1)%N,jSite);
  deltaEnergy -= inplaneEnergy((N+iSite-1)%N,jSite);
  deltaEnergy -= inplaneEnergy(iSite,(jSite+1)%N);
  deltaEnergy -= inplaneEnergy(iSite,(N+jSite-1)%N);

  if(oldCell!=0){
    deltaEnergy -= volumeEnergy(oldCell);
    //deltaEnergy -= perimeterEnergy(oldCell);
    deltaEnergy -= blobularEnergy(oldCell);
  }

  if(newCell!=0){
    deltaEnergy -= volumeEnergy(newCell);
    //deltaEnergy -= perimeterEnergy(newCell);
    deltaEnergy -= blobularEnergy(newCell);
  }

  // Flip it

  lattice[iSite][jSite][0]=newCell;
  adjustVolumes();
  adjustPerimeters();

  // New energy

  deltaEnergy += outplaneEnergy(iSite,jSite);

  deltaEnergy += inplaneEnergy(iSite,jSite);
  deltaEnergy += inplaneEnergy((iSite+1)%N,jSite);
  deltaEnergy += inplaneEnergy((N+iSite-1)%N,jSite);
  deltaEnergy += inplaneEnergy(iSite,(jSite+1)%N);
  deltaEnergy += inplaneEnergy(iSite,(N+jSite-1)%N);

  if(oldCell!=0){
    deltaEnergy += volumeEnergy(oldCell);
    //deltaEnergy += perimeterEnergy(oldCell);
    deltaEnergy += blobularEnergy(oldCell);
  }

  if(newCell!=0){
    deltaEnergy += volumeEnergy(newCell);
    //deltaEnergy += perimeterEnergy(newCell);
    deltaEnergy += blobularEnergy(newCell);
  }

  // Accept or reject the flip

  if( deltaEnergy < 0 ){
    totalEnergy+=deltaEnergy;
    return 1;
  }
  else if( exp(-1.0*beta*deltaEnergy) > (double)rand()/(double)RAND_MAX ){
    totalEnergy+=deltaEnergy;
    return 1;
  }
  else{
    int tmp=newCell;
    newCell=oldCell;
    oldCell=tmp;
    lattice[iSite][jSite][0]=newCell;
    adjustVolumes();
    adjustPerimeters();
    return 0;
  }

}

/*******************************************************************************/
/*** Chooses a spin to flip ***/

void choose()
{
  do{

    int which,thing;

    // Choose a cell perimeter spin to flip

    oldCell = (rand()%numCells)+1;
    int p = rand()%cellPerimeter[oldCell];
    iSite = cellPerimeterList[oldCell][p][0];
    jSite = cellPerimeterList[oldCell][p][1];
    do{
      which = rand()%2;
      thing = 2*(rand()%2)-1;
      if(which==0)
         newCell = lattice[(N+iSite+thing)%N][jSite][0];
      else
        newCell = lattice[iSite][(N+jSite+thing)%N][0];
    }while( oldCell == newCell );

  // Choose to invade or surrender

    if(rand()%2==0){
      if(which==0){
        iSite = (N+iSite+thing)%N;
        int tmp = oldCell;
        oldCell = newCell;
        newCell = tmp;
      }
      else{
        jSite = (N+jSite+thing)%N;
        int tmp = oldCell;
        oldCell = newCell;
        newCell = tmp;
      }
    }

  }while(!maintainsContiguity());

  return;
}

/*******************************************************************************/
/*** Checks if the invasion would cause a cell to break into multiple pieces ***/
/*** Returns FALSE if so (returns TRUE if the flip would maintainContiguity) ***/

bool maintainsContiguity()
{

  int borders[8] = { lattice[(N+iSite-1)%N][(N+jSite-1)%N][0],
                     lattice[(N+iSite-1)%N][jSite][0],
                     lattice[(N+iSite-1)%N][(jSite+1)%N][0],
                     lattice[iSite][(jSite+1)%N][0],
                     lattice[(iSite+1)%N][(jSite+1)%N][0],
                     lattice[(iSite+1)%N][jSite][0],
                     lattice[(iSite+1)%N][(N+jSite-1)%N][0],
                     lattice[iSite][(N+jSite-1)%N][0] };
  /*
	Count how many neighboring spins are in the same cell.
	If there are none, then this is the last spin of that cell.
	We will not let it disappear.
  */

  int totalCellCount = 0;
  for(int n=0; n<8; n++)
    if(borders[n] == oldCell)
      totalCellCount++;   
 
  if(totalCellCount==0)
     return false;

  /*
	The borders array can never be full of cell sites, so by
	starting at a non-cell site, you are guaranteed being able
	to traverse the entire contiguous region without breaks, as
	you aren't starting in the middle of a cell region.

	So, first find the first non-cell spin.
   */

  int index = 0;
  while(borders[index] == oldCell)
    index++;
     
  /*
	Then move along until the site just before the next cell spin.
  */

  while(index < 7 && borders[index+1] != oldCell)
    index++;

  /*
	Starting at the first cell spin, go around until you have
	hit every cell spin.  If you hit a non-cell spin while
	doing this, flipping Site would break contiguity.
	We will not let this happen.
  */

  int inCellCount = 0;
  while (inCellCount < totalCellCount)
  {
    index=(index+1)%8;
    if(borders[index] != oldCell)
      return false;
    inCellCount++;
  }

  return true;
}

/*******************************************************************************/
/*** Adjusts volumes after a spin flip ***/

void adjustVolumes()
{

  if(oldCell!=0){
    cellVolume[oldCell]--;
    int v=0;
    while(cellVolumeList[oldCell][v][0]!=iSite || cellVolumeList[oldCell][v][1]!=jSite)
      v++;
    while(v<cellVolume[oldCell]){
      cellVolumeList[oldCell][v][0]=cellVolumeList[oldCell][v+1][0];
      cellVolumeList[oldCell][v][1]=cellVolumeList[oldCell][v+1][1];
      v++;
    }
  }

  if(newCell!=0){
    cellVolumeList[newCell][cellVolume[newCell]][0]=iSite;
    cellVolumeList[newCell][cellVolume[newCell]][1]=jSite;
    cellVolume[newCell]++;
  }

  return;
}


/*******************************************************************************/
/*** Adjusts perimeters after a spin flip ***/

void adjustPerimeters()
{

  /*
    Site used to be a perimeter spin in oldCell.
    We had better remove it.
  */

  if(oldCell!=0){
    cellPerimeter[oldCell]--;
    int p=0;
    while(cellPerimeterList[oldCell][p][0]!=iSite || cellPerimeterList[oldCell][p][1]!=jSite)
      p++;
    while(p<cellPerimeter[oldCell]){
      cellPerimeterList[oldCell][p][0]=cellPerimeterList[oldCell][p+1][0];
      cellPerimeterList[oldCell][p][1]=cellPerimeterList[oldCell][p+1][1];
      p++;
    }
  }

  /*
    Site is now a perimeter cell in newCell.
    We had better add it.
  */

  if(newCell!=0){
    cellPerimeterList[newCell][cellPerimeter[newCell]][0]=iSite;
    cellPerimeterList[newCell][cellPerimeter[newCell]][1]=jSite;
    cellPerimeter[newCell]++;
  }

  /*
    Now we need to go through each of Site's four neighbors,
    rechecking if they are perimeter cells.
  */

  int iTest,jTest,iTempA,jTempA,iTempB,jTempB,iTempC,jTempC;

  for(int neighbor=0;neighbor<4;neighbor++){

    // tedious setup
    switch(neighbor){

      case 0:
      iTest = (iSite+1)%N;
      jTest = jSite;
      iTempA = (iTest+1)%N;
      jTempA = jTest;
      iTempB = iTest;
      jTempB = (jTest+1)%N;
      iTempC = iTest;
      jTempC = (N+jTest-1)%N;
      break;

      case 1:
      iTest = (N+iSite-1)%N;
      jTest = jSite;
      iTempA = (N+iTest-1)%N;
      jTempA = jTest;
      iTempB = iTest;
      jTempB = (jTest+1)%N;
      iTempC = iTest;
      jTempC = (N+jTest-1)%N;
      break;

      case 2:
      iTest = iSite;
      jTest = (jSite+1)%N;
      iTempA = (iTest+1)%N;
      jTempA = jTest;
      iTempB = (N+iTest-1)%N;
      jTempB = jTest;
      iTempC = iTest;
      jTempC = (jTest+1)%N;
      break;

      case 3:
      iTest  = iSite;
      jTest  = (N+jSite-1)%N;
      iTempA = (iTest+1)%N;
      jTempA = jTest;
      iTempB = (N+iTest-1)%N;
      jTempB = jTest;
      iTempC = iTest;
      jTempC = (N+jTest-1)%N;
      break;

    }

    // only do this for non-air
    if(lattice[iTest][jTest][0]!=0){

      // is surely a perimeter now, but maybe used to be as well
      if(lattice[iTest][jTest][0]==oldCell){

        // if it wasn't before, then add it
        if( lattice[iTempA][jTempA][0]==oldCell &&
            lattice[iTempB][jTempB][0]==oldCell &&
            lattice[iTempC][jTempC][0]==oldCell )
        {
          cellPerimeterList[oldCell][cellPerimeter[oldCell]][0]=iTest;
          cellPerimeterList[oldCell][cellPerimeter[oldCell]][1]=jTest;
          cellPerimeter[oldCell]++;
        }

      }

      // used to be a perimeter, might not be any more
      else if(lattice[iTest][jTest][0]==newCell){

        // if it isn't any more, then remove it
        if( lattice[iTempA][jTempA][0]==newCell &&
            lattice[iTempB][jTempB][0]==newCell &&
            lattice[iTempC][jTempC][0]==newCell )
        {
          cellPerimeter[lattice[iTest][jTest][0]]--;
          int p=0;
          while(cellPerimeterList[newCell][p][0]!=iTest || cellPerimeterList[newCell][p][1]!=jTest)
            p++;
          while(p<cellPerimeter[newCell]){
            cellPerimeterList[newCell][p][0]=cellPerimeterList[newCell][p+1][0];
            cellPerimeterList[newCell][p][1]=cellPerimeterList[newCell][p+1][1];
            p++;
          }
        }

      }

    }

  }

  return;

}

/*******************************************************************************/





























