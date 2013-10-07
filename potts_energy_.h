/*******************************************************************************/
#include <set>
#include <utility>

/*** HAMILTONIAN FUNCTIONS ***/

  double  inplaneEnergy(int,int);
  double  outplaneEnergy(int,int);
  double  Hamiltonian();
  double  interactionEnergy(int);
  double  volumeEnergy(int);
  double  anisotropyEnergy(int);
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
    energy += anisotropyEnergy(cell);
    energy += blobularEnergy(cell);
  }
  return energy;
}

/*******************************************************************************/
/*** Returns the interaction energy of a cell ****/

double interactionEnergy(int cell)
{
  double energy = 0.0;

  for(std::set< std::pair<int, int> >::iterator it = cellPerimeterList[cell].begin(); it!=cellPerimeterList[cell].end(); ++it)
	  energy += inplaneEnergy( it->first , it-> second );

  for(std::set< std::pair<int, int> >::iterator it = cellVolumeList[cell].begin(); it!=cellVolumeList[cell].end(); ++it)
	  energy += outplaneEnergy( it->first , it-> second );

  return energy;
}

/*******************************************************************************/
/*** Returns the volume energy of a cell ***/

double volumeEnergy(int cell)
{
  return L_vol*((double)cellVolumeList[cell].size()-targetVolume)*((double)cellVolumeList[cell].size()-targetVolume);
}

/*******************************************************************************/
/*** Returns the perimeter energy of a cell ***/

double anisotropyEnergy(int cell)
{
		return L_ani*(double)cellPerimeterList[cell].size()/(double)cellVolumeList[cell].size();
}

/*******************************************************************************/
/*** Returns the blobular energy of a cell ***/

double blobularEnergy(int cell)
{

  int energy = 0;
  int number = 0;

  for(std::set< std::pair<int, int> >::const_iterator it1 = cellPerimeterList[cell].begin(); it1!=cellPerimeterList[cell].end(); ++it1){
	  for(std::set< std::pair<int, int> >::const_iterator it2 = cellPerimeterList[cell].begin(); it2!=cellPerimeterList[cell].end(); ++it2){

		  int ai = it1->first;
		  int aj = it1->second;
		  int bi = it2->first;
		  int bj = it2->second;

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

  return L_blb * (double)energy / (double)(cellPerimeterList[cell].size());

}

/*******************************************************************************/


























