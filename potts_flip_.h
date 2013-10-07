#include <utility>
#include <set>
#include <vector>
#include <map>

/*******************************************************************************/


/*** SPIN FLIP FUNCTIONS ***/

  int    flip();
  void   choose();
  void   adjustVolumes(int, int);
  void   adjustPerimeters(int);
  bool   maintainsContiguity();
  std::map< std::pair<int, int> , int > calculateChunkSites(int, int);

  int    iSite;
  int    jSite;
  int    oldCell;
  int    newCell;


/*******************************************************************************/
/*** Flip a spin and accept or reject ***/

int flip()
{

  // Choose a spin
  // The site to be invaded is stored in iSite, jSite
  // newCell invades oldCell

  choose();
  std::map< std::pair<int, int> , int > chunk;
  chunk = calculateChunkSites(iSite, jSite);



  // Subtract out parts of old energy associated with spin site

  double deltaEnergy = 0.0;

  for(std::map< std::pair<int, int> , int >::iterator it = chunk.begin(); it != chunk.end(); ++it){
	  std::pair<int,int> point = it->first;
	  deltaEnergy -= outplaneEnergy( point.first, point.second );
	  deltaEnergy -= inplaneEnergy( point.first, point.second );

	  // the four neighbor directions
	  if( chunk.find( std::make_pair( point.first + 1, point.second ) ) != chunk.end() ){
		  deltaEnergy -= inplaneEnergy( point.first + 1, point.second );
		  deltaEnergy -= outplaneEnergy( point.first + 1, point.second );
	  }

	  if( chunk.find( std::make_pair( (N + point.first - 1)%N, point.second ) ) != chunk.end() ){
		  deltaEnergy -= inplaneEnergy( (N + point.first - 1)%N, point.second );
		  deltaEnergy -= outplaneEnergy( point.first + 1, point.second );
	  }

	  if( chunk.find( std::make_pair( point.first, point.second + 1 ) ) != chunk.end() ){
		  deltaEnergy -= inplaneEnergy( point.first, point.second + 1 );
		  deltaEnergy -= outplaneEnergy( point.first + 1, point.second );
	  }

	  if( chunk.find( std::make_pair( point.first, ( N + point.second - 1)%N ) ) != chunk.end() ){
		  deltaEnergy -= inplaneEnergy( point.first, ( N + point.second - 1)%N );
		  deltaEnergy -= outplaneEnergy( point.first + 1, point.second );
	  }
  }

  if(oldCell!=0){
    deltaEnergy -= volumeEnergy(oldCell);
    deltaEnergy -= anisotropyEnergy(oldCell);
    deltaEnergy -= blobularEnergy(oldCell);
  }

  if(newCell!=0){
    deltaEnergy -= volumeEnergy(newCell);
    deltaEnergy -= anisotropyEnergy(newCell);
    deltaEnergy -= blobularEnergy(newCell);
  }

  // Flip it

  for(std::map< std::pair<int, int> , int >::iterator it = chunk.begin(); it != chunk.end(); ++it){
	  std::pair<int,int> point = it->first;
	  lattice[ point.first ][ point.second ][0]=newCell;
	  adjustVolumes( point.first , point.second );
  }
  adjustPerimeters( newCell );
  adjustPerimeters( oldCell );

  // Add in energy associated with flipped site

  for(std::map< std::pair<int, int> , int >::iterator it = chunk.begin(); it != chunk.end(); ++it){
	  std::pair<int,int> point = it->first;
	  deltaEnergy += outplaneEnergy( point.first, point.second );
	  deltaEnergy += inplaneEnergy( point.first, point.second );

	  // the four neighbor directions
	  if( chunk.find( std::make_pair( point.first + 1, point.second ) ) != chunk.end() )
	  {
		  deltaEnergy += inplaneEnergy( point.first + 1, point.second );
		  deltaEnergy += outplaneEnergy( point.first + 1, point.second );
	  }

	  if( chunk.find( std::make_pair( (N + point.first - 1)%N, point.second ) ) != chunk.end() ){
		  deltaEnergy += inplaneEnergy( (N + point.first - 1)%N, point.second );
		  deltaEnergy += outplaneEnergy( point.first + 1, point.second );
	  }

	  if( chunk.find( std::make_pair( point.first, point.second+1 ) ) != chunk.end() ){
		  deltaEnergy += inplaneEnergy( point.first, point.second+1 );
		  deltaEnergy += outplaneEnergy( point.first + 1, point.second );
	  }

	  if( chunk.find( std::make_pair( point.first, ( N + point.second - 1)%N ) ) != chunk.end() ){
		  deltaEnergy += inplaneEnergy( point.first, ( N + point.second - 1)%N );
		  deltaEnergy += outplaneEnergy( point.first + 1, point.second );
	  }
  }

  if(oldCell!=0){
    deltaEnergy += volumeEnergy(oldCell);
    deltaEnergy += anisotropyEnergy(oldCell);
    deltaEnergy += blobularEnergy(oldCell);
  }

  if(newCell!=0){
    deltaEnergy += volumeEnergy(newCell);
    deltaEnergy += anisotropyEnergy(newCell);
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
  for(std::map< std::pair<int, int> , int >::iterator it = chunk.begin(); it != chunk.end(); ++it){
	  std::pair<int,int> point = it->first;
  	  lattice[ point.first ][ point.second ][0]= it->second;
  	  adjustVolumes( point.first, point.second ); // uses oldCell and newCell
    }
  	  adjustPerimeters( newCell );
  	  adjustPerimeters( oldCell );

    return 0;
  }

}

/*******************************************************************************/
/*** Chooses a spin to flip ***/

void choose()
{

  do{

    int which,thing;

    // Choose a cell
    oldCell = (rand()%numCells)+1;

    // Choose a perimeter site within that cell to either extend or surrender
    std::set< std::pair<int, int> >::const_iterator it(cellPerimeterList[oldCell].begin());
    int p = rand()%cellPerimeterList[oldCell].size();
    advance(it,p);
    iSite = it->first;
    jSite = it->second;


    // Choose a neighboring site either up, down, left, or right
    // newCell invades oldCell
    // bit of a misnomer; either one of these could be air

    do{
      which = rand()%2;
      thing = 2*(rand()%2)-1;
      if(which==0)
         newCell = lattice[(N+iSite+thing)%N][jSite][0];
      else
        newCell = lattice[iSite][(N+jSite+thing)%N][0];
    }while( oldCell == newCell );


	  // Choose to invade or surrender
	  // (that is, if you randomly generate an even number, then
	  // switch from surrendering to invading)

    if(rand()%2==0){
      int tmp = oldCell;
      oldCell = newCell;
      newCell = tmp;
      if(which==0)
        iSite = (N+iSite+thing)%N;
      else
        jSite = (N+jSite+thing)%N;
    }

  }while(!maintainsContiguity());

  return;
}

/*******************************************************************************/
/*** Checks if the invasion would cause a cell to break into multiple pieces ***/
/*** Returns FALSE if so (returns TRUE if the flip would maintainContiguity) ***/

bool maintainsContiguity()
{

  std::vector<int> borders;

  // If the o's in the diagram below are sites in the chunk to the flipped
  // then borders are the x sites

  /* Diagram (chunkSize = 2 here, 0 marks (iSite,jSite) )
  x x x x x x x
  x o o o o o x
  x o o o o o x
  x o o 0 o o x
  x o o o o o x
  x o o o o o x
  x x x x x x x
  */

  /* add the top row of x's
  X X X X X X x
  x o o o o o x
  x o o o o o x
  x o o 0 o o x
  x o o o o o x
  x o o o o o x
  x x x x x x x
  */
  for( int i = iSite - chunkSize - 1; i <= iSite + chunkSize; i++ )
	  borders.push_back( lattice[ (N + i)%N ][ jSite - chunkSize - 1 ][0] );

  /* add the right-side column of x's
  x x x x x x X
  x o o o o o X
  x o o o o o X
  x o o 0 o o X
  x o o o o o X
  x o o o o o X
  x x x x x x x
  */
  for( int j = jSite - chunkSize - 1; j <= jSite + chunkSize; j++ )
	  borders.push_back( lattice[ iSite + chunkSize + 1  ][ (N + j)%N ][0] );

  // the bottom row
  for( int i = iSite + chunkSize + 1; i >= iSite - chunkSize; i-- )
	  borders.push_back( lattice[ (N + i)%N ][ jSite + chunkSize + 1 ][0] );

  // the left-side column
  for( int j = jSite + chunkSize + 1; j >= jSite - chunkSize; j-- )
	  borders.push_back( lattice[ iSite - chunkSize - 1  ][ (N + j)%N ][0] );

  // NOTE: THE BORDER SITES MUST BE ADDED IN CONTINUOUS ORDER FOR THIS
  // ALGORITHM TO WORK.


  /*
  int borders[8] = { lattice[(N+iSite-1)%N][(N+jSite-1)%N][0],
                     lattice[(N+iSite-1)%N][jSite][0],
                     lattice[(N+iSite-1)%N][(jSite+1)%N][0],
                     lattice[iSite][(jSite+1)%N][0],
                     lattice[(iSite+1)%N][(jSite+1)%N][0],
                     lattice[(iSite+1)%N][jSite][0],
                     lattice[(iSite+1)%N][(N+jSite-1)%N][0],
                     lattice[iSite][(N+jSite-1)%N][0] };
  */

  /*
	Count how many neighboring spins are in the same cell.
	If there are none, then this is the last spin of that cell.
	We will not let it disappear.
  */

  int totalCellCount = 0;
  for(int n=0; n<borders.size(); n++)
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
  while(index < borders.size() && borders[index] == oldCell)
    index++;
     
  /*
	Then move along until the site just before the next cell spin.
  */

  while(index < borders.size()-1 && borders[index+1] != oldCell)
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
    index=(index+1)%(borders.size());
    if(borders[index] != oldCell)
      return false;
    inCellCount++;
  }

  return true;
}

/*******************************************************************************/
/*** Find the square chunk of sites around the chosen flip site corresponding to chunkSize ***/

std::map< std::pair<int, int>, int > calculateChunkSites( int i, int j ){
	std::map< std::pair<int, int>, int > chunk;
	for( int x=i-chunkSize; x<=i+chunkSize; x++ ){
		for( int y=j-chunkSize; y<=j+chunkSize; y++){
			std::pair<int, int> p;
			p.first = x;
			p.second = y;

			chunk.insert( std::pair< std::pair<int,int> , int>(p, lattice[ p.first ][ p.second ][0]) );
		}
	}
	return chunk;
}

/*******************************************************************************/
/*** Adjusts volumes after a spin flip ***/

void adjustVolumes(int i, int j)
{

  if(oldCell!=0){
	  std::set< std::pair<int,int> >::iterator it = cellVolumeList[oldCell].find( std::make_pair(i,j) );
	  if( it != cellVolumeList[oldCell].end() )
		cellVolumeList[oldCell].erase( it );
  }

  if(newCell!=0)
    cellVolumeList[newCell].insert( std::make_pair(i,j) );

  return;
}


/*******************************************************************************/
/*** Adjusts perimeters after a spin flip ***/


void adjustPerimeters( int cell ){
	if( cell != 0 ){
	cellPerimeterList[ cell ].clear();
	calculatePerimeter( cell );
	}
}

/*
void adjustPerimeters(int i, int j)
{
yo
  //  Site used to be a perimeter spin in oldCell.
  //  We had better remove it.

   if(oldCell!=0){
	   std::map< std::pair<int,int> >::iterator it = cellPerimeterList[oldCell].find( std::make_pair(i,j) );
	   if( it != cellPerimeterList[oldCell].end() )
	     cellPerimeterList[oldCell].erase( it );
   }


  //  Site is now a perimeter cell in newCell.
  //  We had better add it.

  if(newCell!=0)
    cellPerimeterList[newCell].insert( std::make_pair(i,j) );


  //  Now we need to go through each of Site's four neighbors,
  //  rechecking if they are perimeter cells.


  int iTest,jTest,iTempA,jTempA,iTempB,jTempB,iTempC,jTempC;

  for(int neighbor=0;neighbor<4;neighbor++){

    // tedious setup
    switch(neighbor){

      case 0:
      iTest = (i+1)%N;
      jTest = j;
      iTempA = (iTest+1)%N;
      jTempA = jTest;
      iTempB = iTest;
      jTempB = (jTest+1)%N;
      iTempC = iTest;
      jTempC = (N+jTest-1)%N;
      break;

      case 1:
      iTest = (N+i-1)%N;
      jTest = j;
      iTempA = (N+iTest-1)%N;
      jTempA = jTest;
      iTempB = iTest;
      jTempB = (jTest+1)%N;
      iTempC = iTest;
      jTempC = (N+jTest-1)%N;
      break;

      case 2:
      iTest = i;
      jTest = (j+1)%N;
      iTempA = (iTest+1)%N;
      jTempA = jTest;
      iTempB = (N+iTest-1)%N;
      jTempB = jTest;
      iTempC = iTest;
      jTempC = (jTest+1)%N;
      break;

      case 3:
      iTest  = i;
      jTest  = (N+j-1)%N;
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
          cellPerimeterList[oldCell].insert( std::make_pair(iTest, jTest));
      }

      // used to be a perimeter, might not be any more
      else if(lattice[iTest][jTest][0]==newCell){
        // if it isn't any more, then remove it
        if( lattice[iTempA][jTempA][0]==newCell &&
            lattice[iTempB][jTempB][0]==newCell &&
            lattice[iTempC][jTempC][0]==newCell ){
     	   std::map< std::pair<int,int> >::iterator it = cellPerimeterList[newCell].find( std::make_pair(iTest,jTest) );
     	   if( it != cellPerimeterList[newCell].end() )
     	     cellPerimeterList[newCell].erase( it );
        }
      }

    }

  }

  return;

}
*/
/*******************************************************************************/
