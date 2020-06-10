#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <time.h>
#include<string.h>
#include "p_vec_ten.hh"

using namespace std;
typedef p_vec<> Vec;

#define SQ(x)     ((x) * (x))
const p_ten2 unity(1.,0.,0.,0.,1.,0.,0.,0.,1.);


double incrStr;
int NRPART;
double LX;
double LY;
double LZ;

int t_first;
int num_times;
int t_step;
int t_delta;

double intvl;
double gdot;

const int MAX_NB=200;
const int MAX_NB_Large = 400;
const double CUTOFF_NBLIST = 1.5;		/// should be equal to first min. of pair-correl function
const double CUTOFF_NBLarge = 4.;
const double CUTOFF_NBLIST_SQ = SQ(CUTOFF_NBLIST);
/// subbox parameters
const double BOX_CUTOFF=2.2*CUTOFF_NBLarge;
const int MAX_NRMOL_PER_SUBBOX=3000;

const double LAYER_HWIDTH = 0.5;			/// half-width of layer around z=0 for correl. func. computation

double deltaR_Rdf;		/// resolution
double rangeRdf;
int sizeHistRdf;
int sizeHistAngle;

double LXINV;
double LYINV;
double LZINV;
double XMIN;
double YMIN;
double ZMIN;
double XMAX;
double YMAX;
double ZMAX;


double LX0,LY0,LZ0;
double LXINV0;
double LYINV0;
double LZINV0;
double XMIN0;
double YMIN0;
double ZMIN0;
double XMAX0;
double YMAX0;
double ZMAX0;

int *typen;


struct TSubBox
{
    int nrmol;
    int molptr[MAX_NRMOL_PER_SUBBOX];
};
void periodicDiff(int time, const Vec& p1, const Vec& p2, Vec& dr)
{
  /*dr = p1 - p2;
	if(dr.x > 0.5*LX){
		dr.x -= LX;
	}else if(dr.x < -0.5*LX){
		dr.x += LX;
	}
	if(dr.y > 0.5*LY){
		dr.y -= LY;
	}else if(dr.y < -0.5*LY){
		dr.y += LY;
	}*/
	dr = p1 - p2;
//these part acts as the old periodic difference function
//dr.x -= ((int)(2.*dr.x*LXINV)) * LX ;
dr.z -= ((int)(2.*dr.z*LZINV0)) * LZ0 ;
//dr.z -= ((int)(2.*dr.z*LZINV)) * LZ ;
	//dr.z -= ((int)(2.*dr.z*LZINV)) * LZ ;
}
void apply_pbc(int time, Vec& mol)
{
	mol.z = mol.z - ((int) floor((mol.z - ZMIN) / LZ)) * LZ;
	//mol.z = mol.z - ((int) floor((mol.z - ZMIN) / LZ))* LZ;
	//mol.y = mol.y - ((int) floor((mol.y - YMIN) / LY) ) * LY;
	//if ( mol.x < XMIN ) mol.x += LX ;
  //if ( mol.y < YMIN ) mol.y += LY;
	//if ( mol.z < ZMIN ) mol.z += LZ ;

}
void loadData(int t, Vec* data, Vec* posU)
{
	int i, iniLine= 5, tempInd;
	char fname[1200];
	sprintf(fname, "../notch_%d.dump", t);
	cout << "fname: " << fname << endl;
	ifstream inputFile(fname);
	if (!inputFile)
	{
		cout << "ERROR: no input file: " << fname << endl;
		exit(1);
	}
	string temp;
	for(i = 0; i < iniLine; i++){
		getline(inputFile,temp);
	}
  inputFile >> XMIN >> XMAX ;
  LX = XMAX - XMIN;
  inputFile >> YMIN >> YMAX ;
  LY = YMAX - YMIN;
  inputFile >> ZMIN >> ZMAX ;
  LZ = ZMAX - ZMIN;

  LXINV = 1./LX;
  LYINV = 1./LY;
  LZINV = 1./LZ;

  cout << "Ls: " << LX << " " << LY << " " << LZ << endl;

  getline(inputFile,temp);
  cout << temp<< endl;
  getline(inputFile,temp);
  double xs, ys, zs, tmpD;
  inputFile.precision(15);
	for(i = 0; i < NRPART; i++){
		inputFile >> tempInd;
		inputFile >> typen[tempInd - 1] >>  xs >> ys >> zs >> tmpD >> tmpD >> tmpD;
    if(i == 0) {
      cout << tempInd << " " << typen[tempInd - 1] << " " <<  xs << " " << ys << " " << zs << endl;
    }
    /*if(tempInd <= 0 || tempInd >= NRPART){
      cout << "loading function: the particle index is out of range" << endl;
    }*/
    /*if(xs < 0 || xs > 1){
      cout << "xs: " << xs << endl;
    }
    if(ys < 0 || ys > 1){
      cout << "ys: " << ys << endl;
    }
    if(zs < 0 || zs > 1){
      cout << "zs: " << zs << endl;
    }*/
    data[tempInd - 1].x = xs * LX + XMIN;
    data[tempInd - 1].y = ys * LY + YMIN;
    data[tempInd - 1].z = zs * LZ + ZMIN;
		posU[tempInd - 1] = data[tempInd - 1];
		apply_pbc(t,data[tempInd - 1]);
	}
	inputFile.close();
}

void loadDataFile(char *path,int t, Vec* data, Vec *posU, Vec* vel, int** ind){
	int i, iniLine_1= 22,iniLine_2= 13,iniLine_3 = 4, tempInd;
	char fname[200];
	sprintf(fname, "%s/data.%d",path,t);
	ifstream inputFile(fname);
	if (!inputFile)
	{
		cout << "ERROR: no input file: " << fname << endl;
		exit(1);
	}
	inputFile.precision(25);
	string temp;
	double tempD;
	for(i = 0; i < iniLine_1; i++){
		getline(inputFile,temp);
	}
	/*inputFile >> XMIN0 >> XMAX0 >> temp >> temp;
	inputFile >> YMIN0 >> YMAX0 >> temp >> temp;
	inputFile >> ZMIN0 >> ZMAX0 >> temp >> temp;
	cout << XMIN0 << "\t"<< XMAX0 <<endl;
	cout << YMIN0 << "\t"<< YMAX0 <<endl;
	cout << ZMIN0 << "\t"<< ZMAX0 <<endl;
	for(int i = 0; i < iniLine_2; i++){
		getline(inputFile,temp);
	}
	LX0 = abs(XMIN0 - XMAX0); LY0 = abs(YMIN0 - YMAX0); LZ0 = abs(ZMIN0 - ZMAX0);
	LXINV0 = 1./(double)LX0; LYINV0 = 1./(double)LY0; LZINV0 = 1./(double)LZ0;*/
	for(i = 0; i < NRPART; i++){
		inputFile >> tempInd;
		inputFile >> typen[tempInd - 1] >> data[tempInd - 1].x >> data[tempInd - 1].y >> data[tempInd - 1].z >> ind[tempInd - 1][0]>> ind[tempInd - 1][1]>> ind[tempInd - 1][2];
		posU[tempInd - 1 ].x = data[tempInd - 1].x + ind[tempInd - 1][0] * LX;
		posU[tempInd - 1 ].y = data[tempInd - 1].y + ind[tempInd - 1][1] * LY;
		posU[tempInd - 1 ].z = data[tempInd - 1].z + ind[tempInd - 1][2] * LZ;
		data[tempInd - 1] = posU[tempInd - 1];
		apply_pbc(t, data[tempInd - 1]);
	}
	for(int i = 0; i < iniLine_3; i++){
		getline(inputFile,temp);
		//cout << temp<< endl;
	}
	for(i = 0; i < NRPART; i++){
		inputFile >> tempInd;
		inputFile >> vel[tempInd - 1].x >> vel[tempInd - 1].y >> vel[tempInd - 1].z;

	}
	inputFile.close();
}
void WriteOutData(Vec* pos,int* typenDef, Vec* vel,char* fname){
	//char fname[30];
	//sprintf(fname,"data.mod");
	ofstream oFile(fname);
	oFile.precision(25);
	oFile << "#LAMMPS data file created by muhammad\n" << endl;
	oFile << NRPART << " atoms"<< endl;
	oFile << "8 atom types\n" << endl;
	oFile << XMIN << " " << XMAX << " xlo xhi"<< endl;
	oFile << YMIN << " " << YMAX << " ylo yhi"<< endl;
	oFile << ZMIN << " " << ZMAX << " zlo zhi\n"<< endl;
	//oFile << 0. << " " << 0. << " " << 0. << " xy xz yz\n" << endl;
	oFile << "Masses\n"<< endl;
	oFile << 1 << " " << 1<< endl;
	oFile << 2 << " " << 1<< endl;
	oFile << 3 << " " << 1<< endl;
	oFile << 4 << " " << 1<< endl;
	oFile << 5 << " " << 1<< endl;
	oFile << 6 << " " << 1<< endl;
	oFile << 7 << " " << 1<< endl;
	oFile << 8 << " " << 1<< endl;
	oFile << "\nPair Coeffs\n"<< endl;
	oFile << 1 << " "<< 1 << " "<< 1 <<endl;
	oFile << 2 << " "<< 0.5 << " "<< 0.88 <<endl;
	oFile << 3 << " "<< 1 << " "<< 1 <<endl;
	oFile << 4 << " "<< 0.5 << " "<< 0.88 <<endl;
	oFile << 5 << " "<< 1 << " "<< 1 <<endl;
	oFile << 6 << " "<< 0.5 << " "<< 0.88 <<endl;
	oFile << 7 << " "<<1 << " "<< 1 <<endl;
	oFile << 8 << " "<< 0.5 << " "<< 0.88 <<endl;
	oFile << "\nAtoms\n"<< endl;
	for(int i = 0; i < NRPART; i++){
		oFile << i+1 <<" "<< typenDef[i] << " "<< pos[i] << " 0 0 0" << endl;
	}
	oFile << "\nVelocities\n" << endl;
	for(int i = 0; i < NRPART; i++){
		oFile << i+1 <<" "<< vel[i] << endl;
	}
	oFile.close();
}
double CalPhi(double r){
	double phi;
	phi = 1./(double)pow((double)M_PI,1.5)*exp(-(r*r));
	return(phi);
}
void makeGridBased(int t, Vec *pos, Vec* uField, Vec ***newR, Vec ***uSMonGird, double ***rhoOnGird, double boxShift, int Nx, int Ny, int Nz, double delLoc){

	int ix, iy,iz, i, j, jj, subShift;//, , nrpairs;
	int ix_help, iy_help, iz_help;
	int ixn, iyn, izn, nr_particles_in_the_subbox;

	int NR_SUBX, NR_SUBY, NR_SUBZ;
	double LSUBX, LSUBY, LSUBZ, LSUBXINV, LSUBYINV, LSUBZINV;

	TSubBox ***sub;
	//int period = 1000;
	//ofstream test("testNeighboring.txt");
	/// init subbox structure
    NR_SUBX = (int)ceil(LX0/BOX_CUTOFF);
    LSUBX = LX0/NR_SUBX;
    LSUBXINV = 1.0/LSUBX;

    NR_SUBY = (int)ceil(LY0/BOX_CUTOFF);
    LSUBY = LY0/NR_SUBY;
    LSUBYINV = 1.0/LSUBY;

    NR_SUBZ = (int)ceil(LZ0/BOX_CUTOFF);
    LSUBZ = LZ0/NR_SUBZ;
    LSUBZINV = 1.0/LSUBZ;
    cout << "LSUBs: " <<  LSUBX << " " << LSUBY << " " << LSUBZ << endl;
    cout << "NR_SUBs: "<< NR_SUBX << " " << NR_SUBY << " " << NR_SUBZ << endl;

    subShift = 0.;//(int) floor(boxShift * LSUBXINV);

	  cout << "the sub shift and the box shift: " << subShift << " " << boxShift << endl;


	  sub = (TSubBox***) malloc(NR_SUBX*sizeof(TSubBox**));
    for (i = 0; i < NR_SUBX; i++){
			sub[i] = (TSubBox**) malloc(NR_SUBY*sizeof(TSubBox*));
			for(j = 0; j<NR_SUBY; j++){
				sub[i][j] = (TSubBox*) malloc(NR_SUBZ*sizeof(TSubBox));
			}
    }

	Vec dr;

	for (ix=0; ix<NR_SUBX; ix++)
		for (iy=0; iy<NR_SUBY; iy++)
			for (iz=0; iz<NR_SUBZ; iz++)
				sub[ix][iy][iz].nrmol = 0;

  cout<< "assigning the particles to the boxes" << endl;
	//redistribute particles among subboxes
	for(i=0;  i< NRPART; i++){
		if (pos[i].x < XMIN0 || pos[i].x>XMAX0 || pos[i].y<YMIN0 || pos[i].y>YMAX0 || pos[i].z<ZMIN0 || pos[i].z>ZMAX0)//||
		{
      cout << "In function makeGridBased()" << endl;
			cout << "ERROR: Particle " << i << " is out of simulation box." << endl;
      cout << pos[i] << endl;
      cout << "X ext.: " << XMIN0 << " " << XMAX0 << endl;
      cout << "Y ext.: " << YMIN0 << " " << YMAX0 << endl;
      cout << "Z ext.: " << ZMIN0 << " " << ZMAX0 << endl;
			exit(22);
		}

		ix=(int) floor((pos[i].x - XMIN0 )*LSUBXINV);
		iy=(int) floor((pos[i].y - YMIN0 )*LSUBYINV);
		iz=(int) floor((pos[i].z - ZMIN0 )*LSUBZINV);
		/*if(iz < 0){
			iz = 0;
		}else if(iz >= NR_SUBZ){
			iz = NR_SUBZ - 1;
		}*/
		if(ix < 0 || ix >= NR_SUBX || iy< 0 ||  iy>=NR_SUBY ||  iz < 0 || iz>= NR_SUBZ)
		{
      cout << "In function makeGridBased()" << endl;
			cout << "ERROR: Subbox index out of range." << endl;
			exit(13);
		}
		sub[ix][iy][iz].molptr[sub[ix][iy][iz].nrmol] = i;
		sub[ix][iy][iz].nrmol++;
		//cout << "reached here too" << endl;
		if (sub[ix][iy][iz].nrmol >= MAX_NRMOL_PER_SUBBOX)
		{
			cout << "ERROR: MAX_NRMOL_PER_SUBBOX exceeded.";
			exit(20);
		}

	}
  cout << "creating the latices" << endl;
	for(int i = 0; i < Nx; i++){
		for(int j = 0; j < Ny; j++){
			for(int k = 0; k < Nz; k++){
				newR[i][j][k].x = XMIN0 + (i+0.5)*delLoc;
				newR[i][j][k].y = YMIN0 + (j+0.5)*delLoc;
				newR[i][j][k].z = ZMIN0 + (k+0.5)*delLoc;
			}
		}
	}
	//subShift = (int) floor(boxShift * LSUBXINV);
  cout << "starting with the main loop!" << endl;
  int nr = 0;
	for(int iix = 0; iix < Nx; iix++){
		for(int jjy = 0; jjy < Ny; jjy++){
			for(int kkz = 0; kkz < Nz; kkz++){
				uSMonGird[iix][jjy][kkz].reset();
				rhoOnGird[iix][jjy][kkz] = 0.;
				ix = (int) ((newR[iix][jjy][kkz].x - XMIN0 )*LSUBXINV);
				iy = (int) ((newR[iix][jjy][kkz].y - YMIN0 )*LSUBYINV);
				iz = (int) ((newR[iix][jjy][kkz].z - ZMIN0 )*LSUBZINV);
        if(ix < 0){
            cout << "the index for the grid point and boxing is wrong" << endl;
            cout << "ix: "<< ix << endl;
            exit(10);
        }else if(ix >= NR_SUBX){
          cout << "the index for the grid point and boxing is wrong" << endl;
          cout << "ix: "<< ix << endl;
          exit(10);
        }
        if(iy < 0){
          cout << "the index for the grid point and boxing is wrong" << endl;
          cout << "iy: "<< iy << endl;
          exit(10);
        }else if(iy >= NR_SUBY){
          cout << "the index for the grid point and boxing is wrong" << endl;
          cout << "iy: "<< iy << endl;
          exit(10);
        }
        if(iz < 0){
          cout << "the index for the grid point and boxing is wrong" << endl;
          cout << "iz: "<< iz << endl;
          exit(10);
        }else if(iz >= NR_SUBZ){
          cout << "the index for the grid point and boxing is wrong" << endl;
          cout << "iz: "<< iz << endl;
          exit(10);
        }
				for(ix_help = ix - 1; ix_help <ix + 2; ix_help++){
					for(iy_help = iy - 1; iy_help < iy + 2; iy_help++){
						for(iz_help = iz - 1; iz_help < iz + 2; iz_help++){
              ixn = ix_help;izn = iz_help;iyn =iy_help;
              if(ixn < 0){
                //ixn += NR_SUBX;
                continue;
              }else if(ixn >= NR_SUBX){
                //ixn -= NR_SUBX;
                continue;
              }
              //iyn = iy_help + iy;
              if(iyn < 0){
                //ixn += NR_SUBX;
                continue;
              }else if(iyn >= NR_SUBY){
                //ixn -= NR_SUBX;
                continue;
              }
              //izn = iz_help + iz;
              if(izn < 0){
                izn += NR_SUBZ;
                //continue;
              }else if(izn >= NR_SUBZ){
                izn -= NR_SUBZ;
                //continue;
              }
							nr_particles_in_the_subbox = sub[ixn][iyn][izn].nrmol;
							for(jj = 0; jj < nr_particles_in_the_subbox; jj++){
								j = sub[ixn][iyn][izn].molptr[jj];
                if(j < 0){
                  cout << "error with particle index" << endl;
                  cout << "j: " << j << endl;
                  exit(88);
                }else if( j >= NRPART){
                  cout << "error with particle index" << endl;
                  cout << "j: " << j << endl;
                  exit(88);
                }
								periodicDiff(t,newR[iix][jjy][kkz],pos[j],dr);
								double r = dr.length();
								if(dr.length() < 1.2*delLoc){
									uSMonGird[iix][jjy][kkz] += uField[j]*CalPhi(r);
									rhoOnGird[iix][jjy][kkz] += CalPhi(r);
								}
							}
						}
					}
				}
				if(rhoOnGird[iix][jjy][kkz] != 0){
					uSMonGird[iix][jjy][kkz] = uSMonGird[iix][jjy][kkz] * (1./(double)rhoOnGird[iix][jjy][kkz]);
				}
        nr++;
        /*if(nr % 100000 == 0){
          cout << "finished: " << nr << endl;
        }*/
			}
		}
	}
  cout << "done here" << endl;
	for (i = 0; i < NR_SUBX; i++){
    for (j = 0; j < NR_SUBY; j++){
		    free(sub[i][j]);
    }
    free(sub[i]);
	}
	free(sub);
}
void createCSV(char *outname, Vec ***r, Vec ***u, int lx, int ly, int lz){
  ofstream outfile;
  outfile.open(outname);
  outfile << "X,Y,Z,Ux,Uy,Uz" << endl;
  for(int i = 0; i < lx ; i++){
    for(int j = 0; j < ly; j++){
      for(int k = 0; k < lz; k++){
        outfile << r[i][j][k].x << ","<< r[i][j][k].y << ","<< r[i][j][k].z << ","<< r[i][j][k].z << "," << 0. << "," << 0. << endl;
      }
    }
  }
}
void vtkVecSave(char *outname, Vec ***r, Vec ***u, int lx, int ly, int lz){

  int ii;

  ofstream outfile;

  outfile.open(outname);
  // Write header of the file
  outfile << "# vtk DataFile Version 2.0\n";;
  outfile << "flow_observables\n";
  outfile <<"ASCII\n";
  outfile << "DATASET RECTILINEAR_GRID\n";
  outfile << "DIMENSIONS "<< lx <<" "<< ly<<" "<< lz << endl;
  outfile << "X_COORDINATES " <<  lx << " float\n";

  for(ii = 0; ii < lx; ii++){
      outfile << ii << " ";
  }
  outfile << endl;
  //fprintf(outfile,"\n");
  outfile << "Y_COORDINATES " << ly << " float" << endl;
  //fprintf(outfile,"Y_COORDINATES %d float\n",ly);

  for(ii = 0; ii < ly; ii++){
    outfile << ii << " ";
    //fprintf(outfile,"%d ",ii);
  }
  outfile << endl;
  //fprintf(outfile,"\n");
  outfile << "Z_COORDINATES " << lz << " float\n";
  //fprintf(outfile,"Z_COORDINATES %d float\n",lz);

  for(ii = 0; ii < lz; ii++){
    outfile << ii << " ";
    //fprintf(outfile,"%d ",ii);
  }
  outfile << "\n"<< endl;
  //fprintf(outfile,"\n");
  //fprintf(outfile,"\n");
  outfile << "POINT_DATA " <<lx*ly*lz << endl;
  //fprintf(outfile,"POINT_DATA %d\n",lx*ly*lz);
  outfile << "VECTORS velocity_vector float" << endl;
  //outfile <<"LOOKUP_TABLE default" << endl;
  for(int i = 0; i < lx ; i++ ){
    for(int j = 0; j < ly; j++){
      for(int k = 0; k < lz; k++){
        outfile << u[i][j][k] << endl;
      }
    }
  }
  //fprintf(outfile,"SCALARS density float 1\n");
  //fprintf(outfile,"LOOKUP_TABLE default\n");
  outfile.close();
  //fclose(outfile);
}

void fakeVtkVecSave(char *outname, Vec ***r, Vec ***u, int lx, int ly, int lz){

  int ii;

  ofstream outfile;

  outfile.open(outname);
  // Write header of the file
  outfile << "# vtk DataFile Version 2.0\n";;
  outfile << "flow_observables\n";
  outfile <<"ASCII\n";
  outfile << "DATASET RECTILINEAR_GRID\n";
  outfile << "DIMENSIONS "<< lx <<" "<< ly<<" "<< lz << endl;
  outfile << "X_COORDINATES " <<  lx << " float\n";

  for(ii = 0; ii < lx; ii++){
      outfile << ii << " ";
  }
  outfile << endl;
  //fprintf(outfile,"\n");
  outfile << "Y_COORDINATES " << ly << " float" << endl;
  //fprintf(outfile,"Y_COORDINATES %d float\n",ly);

  for(ii = 0; ii < ly; ii++){
    outfile << ii << " ";
    //fprintf(outfile,"%d ",ii);
  }
  outfile << endl;
  //fprintf(outfile,"\n");
  outfile << "Z_COORDINATES " << lz << " float\n";
  //fprintf(outfile,"Z_COORDINATES %d float\n",lz);

  for(ii = 0; ii < lz; ii++){
    outfile << ii << " ";
    //fprintf(outfile,"%d ",ii);
  }
  outfile << "\n"<< endl;
  //fprintf(outfile,"\n");
  //fprintf(outfile,"\n");
  outfile << "POINT_DATA " <<lx*ly*lz << endl;
  //fprintf(outfile,"POINT_DATA %d\n",lx*ly*lz);
  outfile << "VECTORS velocity_vector float" << endl;
  //outfile <<"LOOKUP_TABLE default" << endl;
  for(int i = 0; i < lx ; i++ ){
    for(int j = 0; j < ly; j++){
      for(int k = 0; k < lz; k++){
        outfile << r[i][j][k].z << " " << 0. << " " << 0. << endl;
      }
    }
  }
  //fprintf(outfile,"SCALARS density float 1\n");
  //fprintf(outfile,"LOOKUP_TABLE default\n");
  outfile.close();
  //fclose(outfile);
}

void StoreUField(int t, char* fname,Vec* pos,Vec *UField){
	ofstream oFile(fname);
	oFile.precision(15);
	oFile << "ITEM: TIMESTEP\n"<< t<<"\nITEM: NUMBER OF ATOMS\n"<<NRPART<<"\nITEM: BOX BOUNDS xy xz yz pp pp pp\n"<< XMIN<< " "<< XMAX <<" 0\n"<< YMIN << " "<< YMAX <<" 0\n"<< ZMIN<< " "<< ZMAX <<" 0"<< endl;
	oFile << "ITEM: ATOMS id type xu yu zu ux uy uz uVal\n";
	for(int i = 0; i < NRPART; i++){
		oFile << i+1 << " "<< typen[i] << " " << pos[i]<< " " << UField[i] << " " << UField[i].length()<< endl;
	}
	oFile.close();
}

void StoreScalarField(int t, char* fname,Vec* pos,double* val){
	ofstream oFile(fname);
	oFile.precision(15);
	oFile << "ITEM: TIMESTEP\n"<< t<<"\nITEM: NUMBER OF ATOMS\n"<<NRPART<<"\nITEM: BOX BOUNDS xy xz yz pp pp pp\n"<< XMIN<< " "<< XMAX <<" 0\n"<< YMIN << " "<< YMAX <<" 0\n"<< ZMIN<< " "<< ZMAX << " 0"<<endl;
																													oFile << "ITEM: ATOMS id type xu yu zu val\n";
	for(int i = 0; i < NRPART; i++){
		oFile << i+1 << " "<< typen[i] << " " << pos[i]<< " " << val[i] << endl;
	}
	oFile.close();
}
void Store2DVec(int t, char* fname,int NGx, int NGz, double delLoc ,Vec** u){
	ofstream oFile(fname);
	oFile.precision(15);
	oFile << "ITEM: TIMESTEP\n"<< t<<"\nITEM: NUMBER OF ATOMS\n"<<NGx*NGz<<"\nITEM: BOX BOUNDS xy xz yz pp pp pp\n"<< XMIN<< " "<< XMAX <<" 0\n"<< -1 << " "<< 1 <<" 0\n"<< ZMIN<< " "<< ZMAX << " 0"<<endl;
	oFile << "ITEM: ATOMS id type x y z ux uz\n";
  int nr = 1;
	for(int i = 0; i < NGx; i++){
    for(int j = 0; j < NGz; j++){
        oFile << nr << " "<< 1 << " " << XMIN+(i+0.5)*delLoc << " 0. "<< ZMIN+(j+0.5)*delLoc << " " << u[i][j].x << " "<<  u[i][j].z << endl;
        nr++;
    }
	}
	oFile.close();
}
void Store2DVec(int t, char* fname,int NGx, int NGy, double delLoc ,Vec** u, double **rho){
	ofstream oFile(fname);
	oFile.precision(15);
	oFile << "ITEM: TIMESTEP\n"<< t<<"\nITEM: NUMBER OF ATOMS\n"<<NGx*NGy<<"\nITEM: BOX BOUNDS xy xz yz pp pp pp\n"<< XMIN<< " "<< XMAX <<" 0\n"<< -1 << " "<< 1 <<" 0\n"<< ZMIN<< " "<< ZMAX << " 0"<<endl;
	oFile << "ITEM: ATOMS id type x y z ux uy uz rho\n";
  int nr = 1;
	for(int i = 0; i < NGx; i++){
    for(int j = 0; j < NGy; j++){
        oFile << nr << " "<< 1 << " " << XMIN+(i+0.5)*delLoc << " "<< YMIN+(j+0.5)*delLoc<< " 0. " << u[i][j].x << " "<<  u[i][j].y << " " << u[i][j].z << " " << rho[i][j]<< endl;
        nr++;
    }
	}
	oFile.close();
}
void Store2DScalar(int t, char* fname,int NGx, int NGz, double delLoc ,double** val){
	ofstream oFile(fname);
	oFile.precision(15);
	oFile << "ITEM: TIMESTEP\n"<< t<<"\nITEM: NUMBER OF ATOMS\n"<<NGx*NGz<<"\nITEM: BOX BOUNDS xy xz yz pp pp pp\n"<< XMIN<< " "<< XMAX <<" 0\n"<< -1 << " "<< 1 <<" 0\n"<< ZMIN<< " "<< ZMAX << " 0"<<endl;
	oFile << "ITEM: ATOMS id type x y z val\n";
  int nr = 1;
	for(int i = 0; i < NGx; i++){
    for(int j = 0; j < NGz; j++){
        oFile << nr << " "<< 1 << " " << XMIN+(i+0.5)*delLoc << " 0. "<< ZMIN+(j+0.5)*delLoc << " " << val[i][j] << endl;
        nr++;
    }
	}
	oFile.close();
}
p_ten2 dyadicProd(Vec r1, Vec r2){
	p_ten2 out;
	out.xx = r1.x*r2.x;out.xy = r1.x*r2.y;out.xz = r1.x*r2.z;
	out.yx = r1.y*r2.x;out.yy = r1.y*r2.y;out.yz = r1.y*r2.z;
	out.zx = r1.z*r2.x;out.zy = r1.z*r2.y;out.zz = r1.z*r2.z;
	return(out);
}
p_ten2 transpos(p_ten2 t1){
	p_ten2 t2 = t1;
	double temp;
	temp = t2.xy; t2.xy = t2.yx; t2.yx = temp;
	temp = t2.xz; t2.xz = t2.zx; t2.zx = temp;
	temp = t2.yz; t2.yz = t2.zy; t2.zy = temp;
	return(t2);
}


Vec GradPhi(Vec dr){
	Vec out;
	double r = dr.length();
	out = dr*(-2./(double)pow((double)M_PI,1.5)*exp(-(r*r)));
	return out;
}
void CalSumGradPhi(int t, Vec *pos, Vec *dPhiI, double cutOff, int** nblist){
	int i , j,p;
	Vec dr;
	for(i = 0; i < NRPART; i++){
		dPhiI[i].reset();
		j = 0;
		while(nblist[i][j] >= 0){
			p = nblist[i][j];
			periodicDiff(t,pos[i],pos[p],dr);
			if(abs(pos[i].z-pos[p].z) < LZ/2.){
				if(dr.length() < cutOff){
						dPhiI[i] += GradPhi(dr);
				}
			}
			j++;
		}
	}
}
void CalSmoothDisp(int t, Vec* pos, Vec* posU, int** nblist, Vec* uList, Vec* uSM, double cutOff, double* rhoLoc ){
	Vec uStar;
	for(int i = 0; i < NRPART; i++){
		uSM[i].reset();
		int j = 0;
		if(posU[i].z >= ZMIN && posU[i].z <= ZMAX){
			uSM[i] += uList[i] * CalPhi(0.);
		}else if(posU[i].z < ZMIN){
			uSM[i].x += -1.*uList[i].x * CalPhi(0.);
			uSM[i].y += uList[i].y * CalPhi(0.);
			uSM[i].z += uList[i].z * CalPhi(0.);
		}else if(posU[i].z > ZMAX ){
			uSM[i].x += -1.*uList[i].x * CalPhi(0.);
			uSM[i].y += uList[i].y * CalPhi(0.);
			uSM[i].z += uList[i].z * CalPhi(0.);
		}


			//it is close to the upper boundary
		j = 0;
		while(nblist[i][j] >= 0){
			int p = nblist[i][j];
			Vec dr;
			periodicDiff(t,pos[i],pos[p],dr);
			double r = dr.length();
			if(r < cutOff){
				if(abs(pos[i].z - pos[p].z)< LZ/2.){
					uSM[i] += uList[p] * CalPhi(r);
				}/*else{
					uSM[i].y += uList[p].y * CalPhi(r);
					uSM[i].z += uList[p].z * CalPhi(r);
					uSM[i].x += -1.*uList[p].x * CalPhi(r);
				}*/
			}
			j++;
		}
		uSM[i] = uSM[i] * (1./rhoLoc[i]);
	}
}
void CalLocDensity(int t, Vec* pos, int** nblist, double* rhoLoc , double cutOff){
		for(int i = 0; i < NRPART; i++){
		rhoLoc[i] = 0.;
		rhoLoc[i] += CalPhi(0.);
		int j = 0;
		while(nblist[i][j] >= 0){
			int p = nblist[i][j];
			Vec dr;
			periodicDiff(t,pos[i],pos[p],dr);
			if(abs(pos[i].z-pos[p].z) < LZ/2.){
				double r = dr.length();
				if(r < cutOff){
					rhoLoc[i] += CalPhi(r);
				}
			}
			j++;
		}
	}
}
void PhiBasedCal(int t, Vec* pos, Vec *posU, Vec* uList, p_ten2* epsField, Vec* uSM, int** nblist, double *rhoLoc, double cutOff){
	char fname[40];
	p_ten2 firstSum;
	p_ten2 *duR;
	duR = (p_ten2*) malloc(NRPART*sizeof(p_ten2));
	cout << "calculating the local densit" << endl;
	CalLocDensity(t, pos, nblist,  rhoLoc, cutOff);
	sprintf(fname, "rhoLoc-%d.dat",t);
	//StoreScalarField(t, fname,pos,rhoLoc);
	cout << "calculating the smooth displacement" << endl;
	CalSmoothDisp(t, pos, posU, nblist, uList, uSM,cutOff, rhoLoc);
	sprintf(fname, "uSM-%d.dat",t);
	StoreUField(t,fname,pos,uSM);

	Vec *dPhiI;
	dPhiI = (Vec*) malloc(NRPART*sizeof(Vec));

	cout << "calculating the sum of phi grad" << endl;
	CalSumGradPhi(t,pos, dPhiI, cutOff, nblist);

	cout << "calculating the strain" << endl;
	Vec dr;
	Vec uStar;
	for(int i = 0; i < NRPART; i++){
		duR[i].reset();
		firstSum.reset();
		epsField[i].reset();
		int j = 0;
		while(nblist[i][j] >= 0){
			int p = nblist[i][j];
			periodicDiff(t,pos[i],pos[p],dr);
			if(dr.length() < cutOff){
				uStar = uList[p];
				if(abs(pos[i].z - pos[p].z)< LZ/2.){
					firstSum = firstSum + dyadicProd(uStar,GradPhi(dr));
				}/*else{
					uStar.x = -1*uList[p].x;
					firstSum = firstSum + dyadicProd(uStar,GradPhi(dr));
				}*/
			}
			j++;
		}
		duR[i] = ( firstSum - dyadicProd(uSM[i],dPhiI[i])) * (1./rhoLoc[i]);
		epsField[i]=(duR[i]+transpos(duR[i]))*0.5;
	}
	free(duR);
	free(dPhiI);
}

void createNeighborlist(int t,Vec* pos, int **nblist, int** nbLarge, double boxShift)
{
	int ix, iy,iz, i, j, jj,  subShift;//, , nrpairs;
	int ix_help, iy_help, iz_help;
	int ixn, iyn, izn, nr_particles_in_the_subbox;

	int NR_SUBX, NR_SUBY, NR_SUBZ;
	double LSUBX, LSUBY, LSUBZ, LSUBXINV, LSUBYINV, LSUBZINV;
	int *nrbonds;
	nrbonds = (int*) malloc(NRPART*sizeof(int));
	TSubBox ***sub;

	//int period = 1000;
	//ofstream test("testNeighboring.txt");
	/// init subbox structure
    NR_SUBX=(int)ceil(LX/BOX_CUTOFF);
    LSUBX=LX/NR_SUBX;
    LSUBXINV=1.0/LSUBX;

    NR_SUBY=(int)ceil(LY/BOX_CUTOFF);
    LSUBY=LY/NR_SUBY;
    LSUBYINV=1.0/LSUBY;

    NR_SUBZ=(int)ceil(LZ/BOX_CUTOFF);
    LSUBZ=LZ/NR_SUBZ;
    LSUBZINV=1.0/LSUBZ;
	subShift = (int) floor(boxShift * LSUBXINV);

	cout << "the sub shift and the box shift: " << subShift << " " << boxShift << endl;
	cout << "the SubBox dimensions: " << NR_SUBX << " "<< NR_SUBY<< " " << NR_SUBZ << endl;

	//sub = new TSubBox**[NR_SUBX];
	sub = (TSubBox***) malloc(NR_SUBX*sizeof(TSubBox**));
    for (i = 0; i < NR_SUBX; i++)
    {
		sub[i] = (TSubBox**) malloc(NR_SUBY*sizeof(TSubBox*));
		for (j = 0; j<NR_SUBY; j++)
		{
			sub[i][j] = (TSubBox*) malloc(NR_SUBZ*sizeof(TSubBox));
		}
    }

	Vec dr;
	for (ix=0; ix<NR_SUBX; ix++)
		for (iy=0; iy<NR_SUBY; iy++)
			for (iz=0; iz<NR_SUBZ; iz++)
				sub[ix][iy][iz].nrmol = 0;

	//redistribute particles among subboxes
	for(i=0;  i< NRPART; i++){
		nrbonds[i] = 0;
		if (pos[i].x < XMIN || pos[i].x>XMAX || pos[i].y<YMIN || pos[i].y>YMAX || pos[i].z < ZMIN || pos[i].z > ZMAX)
		{
			cout << "ERROR: Particle " << i << " is out of simulation box." << endl;
			exit(22);
		}

		ix=(int) floor((pos[i].x - XMIN )*LSUBXINV);
		iy=(int) floor((pos[i].y - YMIN )*LSUBYINV);
		iz=(int) floor((pos[i].z - ZMIN )*LSUBZINV);
		/*if(iz < 0){
			iz = 0;
		}else if(iz >= NR_SUBZ){
			iz = NR_SUBZ - 1;
		}*/
		if(ix<0 || ix>NR_SUBX || iy<0 ||  iy>NR_SUBY ||  iz<0 || iz>NR_SUBZ)
		{
			cout << "ERROR: Subbox index out of range." << endl;
			cout << "i: "  << i << " ,pos:" << pos[i] << endl;
			cout << "indices: " << ix << " " << iy << " " << iz << endl;
			exit(13);
		}else if(ix == NR_SUBX){
			ix--;
		}else if(iz == NR_SUBZ){
			iz--;
		}else if(iy == NR_SUBY){
			iy--;
		}
		//cout << "reached here" << endl;
		sub[ix][iy][iz].molptr[sub[ix][iy][iz].nrmol] = i;
		sub[ix][iy][iz].nrmol++;
		//cout << "reached here too" << endl;
		if (sub[ix][iy][iz].nrmol >= MAX_NRMOL_PER_SUBBOX)
		{
			cout << "ERROR: MAX_NRMOL_PER_SUBBOX exceeded.";
			exit(20);
		}

	}
	cout << "the first loop got finished" << endl;
	vector<int> nbPerPart(NRPART),nbPerPartL(NRPART);
	//ofstream oFile("NrNeigh.dat");

	for(i = 0; i < NRPART; i++){
		//nr_nb = 0;
		nbPerPart[i] = 0;
		nbPerPartL[i]=0;
		for(j = 0; j < MAX_NB; j++){
			nblist[i][j] = -1;
		}
		for(j = 0; j < MAX_NB_Large; j++){
			nbLarge[i][j]= -1;
		}
		ix=(int) floor((pos[i].x - XMIN )*LSUBXINV);
		iy=(int) floor((pos[i].y - YMIN )*LSUBYINV);
		iz=(int) floor((pos[i].z - ZMIN )*LSUBZINV);
		if(ix == NR_SUBX){
			ix--;
		}else if(iz == NR_SUBZ){
			iz--;
		}else if(iy == NR_SUBY){
			iy--;
		}
		for (ix_help=ix-1; ix_help<ix+2; ix_help++){
			for (iy_help = 0; iy_help < 2; iy_help++){
				for (iz_help=iz-1; iz_help<iz+2; iz_help++){

					ixn = ix_help;izn = iz_help;iyn =iy_help;
					if(ixn < 0){
						ixn += NR_SUBX;
					}else if(ixn >= NR_SUBX){
						ixn -= NR_SUBX;
					}
					iyn = iy_help + iy;
					if(iy == 1 && iy_help == 1){
						iyn = 0;
					}
					if(izn < 0){
						//izn += NR_SUBZ;
						continue;
					}else if(izn >= NR_SUBZ){
						//izn -= NR_SUBZ;
						continue;
					}

					nr_particles_in_the_subbox=sub[ixn][iyn][izn].nrmol;
					for(jj = 0; jj < nr_particles_in_the_subbox; jj++){
						j =sub[ixn][iyn][izn].molptr[jj];
						if(i != j){
							periodicDiff(0,pos[i],pos[j],dr);
							//dr = PerDiff(pos[i],pos[j],boxShift);
							if(dr.length() <= CUTOFF_NBLIST){
								if(nbPerPart[i] >= MAX_NB){
									cout << "Number of neighbors exceeds the assumed number" << endl;
									exit(10);
								}
								nblist[i][nbPerPart[i]]= j;
								nbPerPart[i]++;
							}
							if(dr.length() <= CUTOFF_NBLarge){
								if(nbPerPartL[i] >= MAX_NB_Large){
									cout << "Number of neighbors exceeds the assumed number; larger" << endl;
									exit(10);
								}
								nbLarge[i][nbPerPartL[i]]= j;
								nbPerPartL[i]++;
							}
						}
					}
				}
			}
		}
		if(i % 100000 == 0){
			cout << "the list built for " << i << endl;
		}
	}
	cout << "freeing some local memories" << endl;
	for (i = 0; i < NR_SUBX; i++)
  {
		for (j = 0; j < NR_SUBY; j++)
		{
			free(sub[i][j]);
		}
		free(sub[i]);
  }
	free(sub);free(nrbonds);
}
void StoreStrain(int t, char* fname, Vec* pos,p_ten2* strain){
	ofstream oFile(fname);
	oFile.precision(15);
	oFile << "ITEM: TIMESTEP\n"<< t<<"\nITEM: NUMBER OF ATOMS\n"<<NRPART<<"\nITEM: BOX BOUNDS xy xz yz pp pp pp\n"<< XMIN<< " "<< XMAX <<" 0\n"<< YMIN << " "<< YMAX <<" 0\n"<< ZMIN<< " "<< ZMAX <<" 0"<<endl;
	oFile << "ITEM: ATOMS id type xu yu zu ";
	for(int i = 0; i < 9; i++){
		oFile << "val[" << i+1<<"] ";
		if(i == 8){
			oFile << endl;
		}
	}
	for(int i = 0; i < NRPART; i++){
		oFile << i+1 << " "<< typen[i] << " " << pos[i] << " "<< strain[i] << endl;
	}
	oFile.close();
}
void storeGBVec(char *fname, int t, Vec ***r, Vec ***u, int Nx, int Ny, int Nz){
	ofstream oFile(fname);
	oFile.precision(15);
	oFile << "X,Y,Z,Ux,Uy,Uz\n";
	//int nr = 0;
	for(int i = 0; i < Nx; i++){
		for(int j = 0; j < Ny; j++){
			for(int k = 0; k < Nz; k++){
					oFile << r[i][j][k].x << ","<< r[i][j][k].y << ","<< r[i][j][k].z << ","<< u[i][j][k].x <<",0.," << u[i][j][k].z << endl;//" " << UField[i].length()

			}
		}

	}
	oFile.close();
}
void storeGBVecForOvito(char *fname, int t, Vec ***r, Vec ***u, int Nx, int Ny, int Nz){
  ofstream oFile(fname);
	oFile.precision(15);
	oFile << "ITEM: TIMESTEP\n"<< t <<"\nITEM: NUMBER OF ATOMS\n"<<Nx*Ny*Nz<<"\nITEM: BOX BOUNDS xy xz yz pp pp pp\n"<< XMIN0<< " "<< XMAX0 <<" 0\n"<< YMIN0 << " "<< YMAX0 <<" 0\n"<< ZMIN0<< " "<< ZMAX0 <<" 0"<< endl;
	oFile << "ITEM: ATOMS id type x y z ux uy uz\n";
	int nr = 0;
	for(int i = 0; i < Nx; i++){
		for(int j = 0; j < Ny; j++){
			for(int k = 0; k < Nz; k++){
					oFile << nr << " 1 " << r[i][j][k].x << " "<< r[i][j][k].y << " "<< r[i][j][k].z << " "<< u[i][j][k].x<< " " << u[i][j][k].y<< " " << u[i][j][k].z << endl;//" " << UField[i].length()
          nr++;
			}
		}

	}
	oFile.close();
}
void storeGBVecForOvito(char *fname, int t, Vec ***r, Vec ***u, double ***rho, int Nx, int Ny, int Nz){
  ofstream oFile(fname);
	oFile.precision(15);
	oFile << "ITEM: TIMESTEP\n"<< t <<"\nITEM: NUMBER OF ATOMS\n"<<Nx*Ny*Nz<<"\nITEM: BOX BOUNDS xy xz yz pp pp pp\n"<< XMIN0<< " "<< XMAX0 <<" 0\n"<< YMIN0 << " "<< YMAX0 <<" 0\n"<< ZMIN0<< " "<< ZMAX0 <<" 0"<< endl;
	oFile << "ITEM: ATOMS id type x y z ux uy uz rho\n";
	int nr = 0;
	for(int i = 0; i < Nx; i++){
		for(int j = 0; j < Ny; j++){
			for(int k = 0; k < Nz; k++){
					oFile << nr << " 1 " << r[i][j][k].x << " "<< r[i][j][k].y << " "<< r[i][j][k].z << " "<< u[i][j][k].x<< " " << u[i][j][k].y<< " " << u[i][j][k].z << " "<<  rho[i][j][k] << endl;//" " << UField[i].length()
          nr++;
			}
		}

	}
	oFile.close();
}
void Project2D(char *fname, int t, Vec ***r, Vec ***u, int Nx, int Ny, int Nz){
  ofstream oFile(fname);
	oFile.precision(15);
	oFile << "X,Z,Ux,Uz\n";
  double valx, valz;
	for(int i = 0; i < Nx; i++){
			for(int k = 0; k < Nz; k++){
        valx = 0.;
        valz = 0.;
        for(int j = 0; j < Ny; j++){
          valx += u[i][j][k].x;
          valz += u[i][j][k].z;
        }
        oFile << r[i][0][k].x << ","<< r[i][0][k].z << ","<< valx/(double)Ny << "," << valz/(double)Ny << endl;
		}
	}
	oFile.close();
}
void Make2DGrid(Vec ***u, Vec **v,double ***rho3D ,double **rho2D, int Nx, int Ny, int Nz){
  for(int i = 0; i < Nx; i++){
			for(int j = 0; j < Ny; j++){
        v[i][j].x = 0.;
        v[i][j].z = 0.;
        v[i][j].y = 0.;
        rho2D[i][j] = 0.;
        for(int k = 0; k < Nz; k++){
          v[i][j].x += u[i][j][k].x;
          v[i][j].z += u[i][j][k].z;
          v[i][j].y += u[i][j][k].y;
          rho2D[i][j] += rho3D[i][j][k];
        }
        v[i][j].x = v[i][j].x / (double)Nz;
        v[i][j].z = v[i][j].z / (double)Nz;
        v[i][j].y = v[i][j].y / (double)Nz;
        rho2D[i][j] = rho2D[i][j] / (double)Nz;
		}
	}
}
void CalGradUyx(Vec **v, double **dUyx, int Nx, int Ny){
  double delLoc = 1.;
  for(int i = 0; i < Nx; i++){
    int iMinus = i-1;
    int iPlus = i+1;
    if(iMinus == -1){
      iMinus = 0;
      for(int j = 0; j < Ny; j++){
        dUyx[i][j] = 1./(delLoc)*(v[iPlus][j].y - v[iMinus][j].y);
      }
    }else if(iPlus == Nx){
      iPlus = Nx-1;
      for(int j = 0; j < Ny; j++){
        dUyx[i][j] = 1./(delLoc)*(v[iPlus][j].y - v[iMinus][j].y);
      }
    }else{
      for(int j = 0; j < Ny; j++){
        dUyx[i][j] = 0.5/(delLoc)*(v[iPlus][j].y - v[iMinus][j].y);
      }
    }
  }
}
void CalGradUxy(Vec **v, double **dUxy, int Nx, int Ny){
  double delLoc = 1.;
  for(int j = 0; j < Ny; j++){
    int jMinus = j-1;
    int jPlus = j+1;

    if(jMinus == -1){
      for(int i = 0; i < Nx; i++){
        dUxy[i][j] = 1./(delLoc)*(v[i][jPlus].x - v[i][j].x);
      }
    }else if(jPlus == Ny){
      for(int i = 0; i < Nx; i++){
        dUxy[i][j] = 1./(delLoc)*(v[i][j].x - v[i][jMinus].x);
      }
    }else{
      for(int i = 0; i < Nx; i++){
        dUxy[i][j] = 0.5/(delLoc)*(v[i][jPlus].x - v[i][jMinus].x);
      }
    }
  }
}
void calShearStr2D(Vec **v, double **dUxy, double **dUyx, double **eps,int Nx, int Ny){
  cout << "calculating yx component" << endl;
  CalGradUyx(v,dUyx, Nx,Ny);
  cout << "calculating xy component" << endl;
  CalGradUxy(v,dUxy, Nx,Ny);
  cout << "summing up" << endl;
  for(int i = 0; i < Nx; i++){
    for(int j = 0; j < Ny; j++){
      eps[i][j] = 0.5 * (dUyx[i][j]+ dUxy[i][j]);
    }
  }
}
void store2DShearComponents(int t , char *fname, double **dUxy, double **dUyx, double **eps, double **rho, int Nx, int Ny, double delLoc ){
  ofstream oFile(fname);
	oFile.precision(15);
	oFile << "ITEM: TIMESTEP\n"<< t<<"\nITEM: NUMBER OF ATOMS\n"<<Nx*Ny<<"\nITEM: BOX BOUNDS xy xz yz pp pp pp\n"<< XMIN<< " "<< XMAX <<" 0\n"<< -1 << " "<< 1 <<" 0\n"<< ZMIN<< " "<< ZMAX << " 0"<<endl;
	oFile << "ITEM: ATOMS id type x y z Uxy Uyx Exy rho\n";
  int nr = 1;
	for(int i = 0; i < Nx; i++){
    for(int j = 0; j < Ny; j++){
        oFile << nr << " "<< 1 << " " << XMIN+(i+0.5)*delLoc << " "<< YMIN+(j+0.5)*delLoc<< " 0. " << dUxy[i][j]<< " "<<  dUyx[i][j] << " " << eps[i][j] << " " << rho[i][j]<< endl;
        nr++;
    }
	}
	oFile.close();

}
void IdentifySB(double **dUzx,int **actFlag, int Nx, int Nz){
  for(int i = 0; i < Nx; i++){
    for(int k = 0; k < Nz; k++){
      if(abs(dUzx[i][k]) > 0.01){
        actFlag[i][k]=1;
      }else{
        actFlag[i][k]=0;
      }
    }
  }
}
void calCorr(double **val,int **flag, double *corr, int *nrCorr, int Nx, int Nz){
  for(int i = 0; i < (int)ceil(Nx/2.); i++){
    nrCorr[i]=0;
    corr[i]=0.;
  }
  double *val1D;
  val1D = (double*) malloc(Nx*sizeof(double));
  int nr = 0;
  for(int i=0 ; i < Nx; i++){
    val1D[i]=0.;
    nr=0;
    for(int k = 0; k < Nz; k++){
      if(flag[i][k] == 1){
        val1D[i] += val[i][k];
        nr++;
      }
    }
    if(nr != 0){
      val1D[i]= val1D[i]/(double)nr;
    }
  }
  int dx;
  for (int i = 0; i < Nx; i++){
    for(int j = i; j < Nx; j++){
      dx= j - i;
      if(dx > int(Nx/2.)){
        dx -= (int)Nx/2.;
      }else if(dx < 0){
        cout << "wrong assumption" << endl;
        exit(859);
      }
      corr[dx] += val1D[i]*val1D[j];
      nrCorr[dx]++;
    }
  }
  for(int i = 0; i < (int)ceil(Nx/2.); i++){
    if(nrCorr[i]!=0){
      corr[i]=corr[i]/(double)nrCorr[i];
      nrCorr[i]=1;
    }
  }
  free(val1D);
}
void Save1DCorr(char *fname, double *corr, int *nrCorr,int sizeN, double delLoc){
  ofstream oFile(fname);
  for(int i = 0; i < sizeN; i++){
    if(nrCorr[i] != 0){
      corr[i] = corr[i] / (double)nrCorr[i];
    }
    oFile << (i)*delLoc << " " << corr[i] << endl;
  }
  oFile.close();
}
void rmAffinePart(Vec **uNon,Vec **v, int Nx, int Nz){
  int nr = 0;
  for(int i = 0; i < Nx; i++){
    for(int k = 0; k < Nz; k++){
      uNon[i][k] = v[i][k];
      int nrFor = 0, nrRew = 0.;
      for( int n = i+1; n < Nx; n++){
        if(v[i][k].x * v[n][k].x >= 0){
          nrFor++;
        }else{
          break;
        }
      }
      for( int n = i-1; n >= 0; n--){
        if(v[i][k].x*v[n][k].x >= 0){
          nrRew++;
        }else{
          break;
        }
      }
      if((nrFor+nrRew) > 125){
        double val = 0.;
        for(int n = i - 1; n >= 0; n--){
          val += v[i][k].x;
        }
        for(int n = 0; n < Nx; n++){
          val += v[i][k].x;
        }
        val /= (double)(nrRew + nrFor + 1.);
        uNon[i][k].x = uNon[i][k].x - val;
      }
      if(nr % 10000 == 0){
        cout << "removing done for: " << nr << endl;
      }
      nr++;
    }
  }
}
void Save2D_disp(char *fname, Vec ***r, Vec** v, int Nx, int Nz){
  ofstream oFile(fname);
	oFile.precision(15);
	oFile << "X,Z,Ux,Uz\n";
  //double valx, valz;
	for(int i = 0; i < Nx; i++){
			for(int k = 0; k < Nz; k++){
        oFile << r[i][0][k].x << ","<< r[i][0][k].z << ","<< v[i][k].x << "," << v[i][k].z << endl;
		}
	}
	oFile.close();
}
void storeThirdGBVec(char *fname, int t, Vec ***r, Vec ***u, int Nx, int Ny, int Nz){
  char fileName[1000];
  for(int f = 0; f < 3; f++){
    sprintf(fileName, "%s-%fpart-%d.csv",fname,f,t);
    ofstream oFile(fileName);
  	oFile.precision(15);
  	oFile << "X,Y,Z,Ux,Uy,Uz\n";
  	//int nr = 0;
  	for(int i = int((f)*Nx/3.); i < min(int((f+1)*Nx/3.),Nx); i++){
  		for(int j = 0; j < Ny; j++){
  			for(int k = 0; k < Nz; k++){
  					oFile << r[i][j][k].x << ","<< r[i][j][k].y << ","<< r[i][j][k].z << ","<< u[i][j][k].x <<",0.," << u[i][j][k].z << endl;//" " << UField[i].length()
  			}
  		}
    }
    oFile.close();
	}

}
p_ten2 computeMeanStr(p_ten2* str){
  p_ten2 out;
  out.reset();
  for(int i = 0; i < NRPART; i++){
    out = out + str[i];
  }
  out = out * (1./(double)NRPART);
  return(out);
}
double computeMean(double* val){
  double out;
  for(int i = 0; i < NRPART; i++){
    out += val[i];
  }
  out = out / (double)NRPART;
  return(out);
}
void CalcFalkStr(int t ,Vec* pos1, Vec* pos2, double *listD2, p_ten2* epsList,int** nblist){
	for(int i = 0; i < NRPART; i++){
		listD2[i] = ((double)rand()/(double)RAND_MAX)*2-1;
	}
	//char fname[100];
	//sprintf(fname,"D2-%d.dat",t);
	//StoreScalarField(t, fname,pos1,listD2);
	//char fname[100];
	int i , j;
	double D2;
	p_ten2 X,Y, epsilon;
	Vec dist_then,dist_now, tmpV11, tmpV12,  tmpV21, tmpV22;
	for (i=0; i < NRPART; i++){
		X.reset();
		Y.reset();
		tmpV11 = pos1[i];
		tmpV12 = pos2[i];
		epsilon.reset();
		j=0;
		while (nblist[i][j] >=0)
		{
			int  p = nblist[i][j];
			tmpV21 = pos1[p];tmpV22=pos2[p];
			periodicDiff(t,tmpV11,tmpV21, dist_then);
			periodicDiff(t, tmpV12, tmpV22,dist_now);
			X.xx += dist_now.x * dist_then.x;
			X.xy += dist_now.x * dist_then.y;
			X.xz += dist_now.x * dist_then.z;
			X.yx += dist_now.y * dist_then.x;
			X.yy += dist_now.y * dist_then.y;
			X.yz += dist_now.y * dist_then.z;
			X.zx += dist_now.z * dist_then.x;
			X.zy += dist_now.z * dist_then.y;
			X.zz += dist_now.z * dist_then.z;
			Y.xx += dist_then.x * dist_then.x;
			Y.xy += dist_then.x * dist_then.y;
			Y.xz += dist_then.x * dist_then.z;
			Y.yx += dist_then.y * dist_then.x;
			Y.yy += dist_then.y * dist_then.y;
			Y.yz += dist_then.y * dist_then.z;
			Y.zx += dist_then.z * dist_then.x;
			Y.zy += dist_then.z * dist_then.y;
			Y.zz += dist_then.z * dist_then.z;
			j++;
		}
		epsilon = X * (Y.invert()) - unity;
		epsList[i] = epsilon;
		D2 = 0;
		j = 0;
		while (nblist[i][j] >= 0)
		{
			int p=nblist[i][j];
			periodicDiff( t, pos1[i], pos1[p], dist_then);
			periodicDiff( t+ t_delta, pos2[i], pos2[p], dist_now);
			D2 += SQ(dist_now - (epsilon + unity) * dist_then);
			j++;
		}
		if(j != 0){
			listD2[i] = D2/(double)j;
		}
	}
	cout << "mean Falk strain: " << computeMeanStr(epsList) << endl;
	cout << "mean D2: " << computeMean(listD2) << endl;
	//cout << "storing the D2 list" << endl;
	//sprintf(fname,"D2-%d.dat",t);
	//StoreScalarField(t, fname,pos1,listD2);
	//free(epsList);
}
int MakeDUzx1D(int **actFlag, double **dUzx, double *dUzx1D,int *nr_dUzx1D, int sizeX, int sizeZ){
  for(int i = 0; i < sizeX; i++){
    dUzx1D[i] = 0.;
    nr_dUzx1D[i] = 0;
  }
  for(int i = 0; i < sizeX; i++){
    for(int j = 0; j < sizeZ; j++){
      if(actFlag[i][j] == 1){
        int ix = i;
        dUzx1D[ix] += dUzx[i][j];
        nr_dUzx1D[ix]++;
      }
    }
  }
  for(int i = 0; i < sizeX; i++){
    if(nr_dUzx1D[i] != 0){
      dUzx1D[i] = dUzx1D[i] / nr_dUzx1D[i];
    }
  }
}


void calcHistDis(int t, double *dUzx1D, int *nr_dUzx1D, double *hist, int sizeX, int sizeZ, double delLoc, int sizeHistDis){
  double *dist;
  dist = (double*) malloc(sizeX * sizeof(double));
  int nrDis = 0;
  //int breakFlag=0;
  int start_flag = 0;
  int x1;
  for(int i = 0; i < sizeX; i++){
    if(start_flag == 0){
      if(nr_dUzx1D[i] != 0){
        x1 = i;
        start_flag = 1;
      }
    }else{
      if(nr_dUzx1D[i] != 0){
        if(dUzx1D[i] * dUzx1D[x1] <= 0){
          dist[nrDis] = (i - x1 + 1) * delLoc;
          nrDis++;
          x1 = i;
        }
      }
    }
  }

  cout << "nrDis: " << nrDis << endl;
  ofstream test("testDis.dat");
  for(int i = 0;  i < nrDis; i++){
    test << i << " " << dist[i] << endl;
  }
  test.close();

  double minDis = 0, maxDis = 200;
  double delHist = 1.;
  for(int i = 0; i < nrDis; i++){
    if(dist[i] >= maxDis || dist[i] < minDis){
        cout << "the estimate for min and max dist. is not correct" << endl;
        cout << i << " " << dist[i] << endl;
        exit(1212);
    }
    int bin = (int)floor((dist[i] - minDis)/delHist);
    if(bin < 0 || bin >= sizeHistDis){
      cout << "error in binning " << endl;
      cout << "bin: " << bin << endl;
      cout << dist[i] << endl;
      exit(222);
    }
    hist[bin]++;
  }
  char fname[100];
  sprintf(fname,"Hist-%d.dat",t);
  test.open(fname);
  for(int i = 0; i < sizeHistDis; i++){
    test << minDis+(i+0.5)*delHist << " " << hist[i] << endl;
  }
  test.close();
  free(dist);
}

int main(int argc, char **argv)
{
	clock_t start, end;
	start = clock();
  if(argc !=5){
    	cout << "the inputs are not correct; their number does not match to the expectations" << endl;
      return 1;

  }

  //char path[1000];
  //LX = atof(argv[1]);
  //LY = atof(argv[2]);
  //LZ = atof(argv[3]);
  NRPART = atoi(argv[1]);
  t_first = atoi(argv[2]);
  t_step = atoi(argv[3]);
  num_times = atoi(argv[4]);
  //intvl = atof(argv[4]);

  cout << "NRPART: "<< NRPART << endl;

   typen = (int*) malloc(NRPART*sizeof(int));

   char fname[100];
   /*int **nblist, **Largenblist;
   nblist = (int**) malloc(NRPART*sizeof(int*));
   Largenblist = (int**) malloc(NRPART*sizeof(int*));
   for(int i = 0; i < NRPART; i++){
  	  nblist[i] = (int*) malloc(MAX_NB*sizeof(int));
  	  Largenblist[i] = (int*) malloc(MAX_NB_Large*sizeof(int));
    }*/

    Vec *pos,*posU, *pos2,*posU2,*uTot;
    int **boxInd, **boxInd2;
    double *listD2;


    listD2 = (double*) malloc(NRPART*sizeof(double));

    p_ten2* epsList;
  	epsList = (p_ten2*) malloc(NRPART*sizeof(p_ten2));

    uTot = (Vec*) malloc(NRPART*sizeof(Vec));

    pos = (Vec*) malloc(NRPART*sizeof(Vec));
    posU = (Vec*) malloc(NRPART*sizeof(Vec));

    pos2 = (Vec*) malloc(NRPART*sizeof(Vec));
    posU2 = (Vec*) malloc(NRPART*sizeof(Vec));

		int NGx,NGy,NGz;
    double delLoc = 3.0;

		Vec ***uSMonGird, ***newR;
		double ***rhoOnGird, **rho2D;
    Vec **v2D;


    int **actFlag;
    double **dUyx,**dUxy, **eps2D;


    double *dUzx1D;
    int *nr_dUzx1D;



    double *corrDu_sum;
    int *nrCorrDu_sum;


    double *corrDu;
    int *nrCorrDu;



    Vec **uNon;


    /*double minDis = 0, maxDis = 200;
    double delHistDis = 1.;
    int sizeHistDis = (int) ceil((maxDis - minDis)/delHistDis);
    double *histDis;
    histDis = (double*) malloc(sizeHistDis * sizeof(double));
    for(int i = 0;  i < sizeHistDis ; i++){
      histDis[i] = 0.;
    }*/


    for(int t = 0; t < num_times; t++){
      int tCurrent = t_first+t*t_step;
      cout << "loading data @ " << tCurrent << endl;
      loadData(tCurrent, pos, posU);
      LX0=LX;LY0=LY;LZ0=LZ;
      LXINV0=LXINV;LYINV0=LYINV;LZINV0=LZINV;
      XMIN0=XMIN;YMIN0=YMIN;ZMIN0=ZMIN;
      XMAX0=XMAX;YMAX0=YMAX;ZMAX0=ZMAX;
      NGx = (int) floor(LX /(double)delLoc);
  		NGy = (int) floor(LY /(double)delLoc);
  		NGz = (int) floor(LZ /(double)delLoc);

      cout << "NX: " << NGx << endl;
      cout << "NY: " << NGy << endl;
      cout << "NZ: " << NGz << endl;

      cout << "the total number: " << NGx*NGy*NGz << endl;

      cout << "allocating the arrays" << endl;
      dUyx = (double**)malloc(NGx*sizeof(double*));
      dUxy = (double**)malloc(NGx*sizeof(double*));
      eps2D = (double**)malloc(NGx*sizeof(double*));
      actFlag = (int**)malloc(NGx*sizeof(int*));

      for(int i = 0; i < NGx; i++){
        dUyx[i] = (double*)malloc(NGy*sizeof(double));
        dUxy[i] = (double*)malloc(NGy*sizeof(double));
        eps2D[i] = (double*)malloc(NGy*sizeof(double));
        actFlag[i] = (int*)malloc(NGy*sizeof(int));
      }



      dUzx1D = (double*) malloc(NGx*sizeof(double));
      nr_dUzx1D = (int*) malloc(NGx*sizeof(double));



      corrDu_sum = (double*) malloc(((int)ceil(NGx/2.))*sizeof(double));
      nrCorrDu_sum = (int*) malloc(((int)ceil(NGx/2.))*sizeof(int));



      for(int i = 0; i < (int)ceil(NGx/2.); i++){
        corrDu_sum[i]=0.;
        nrCorrDu_sum[i]=0.;
      }

      corrDu = (double*) malloc(((int)ceil(NGx/2.))*sizeof(double));
      nrCorrDu = (int*) malloc(((int)ceil(NGx/2.))*sizeof(int));



      uSMonGird= (Vec***) malloc(NGx*sizeof(Vec**));
  		newR = (Vec***) malloc(NGx*sizeof(Vec**));
  		rhoOnGird = (double***) malloc(NGx*sizeof(double**));

      uNon = (Vec**) malloc(NGx*sizeof(Vec*));
      v2D = (Vec**) malloc(NGx*sizeof(Vec*));
      rho2D = (double**) malloc(NGx*sizeof(double*));

  		for(int i = 0; i < NGx; i++){
        v2D[i] = (Vec*) malloc(NGy*sizeof(Vec));
        uNon[i] = (Vec*) malloc(NGy*sizeof(Vec));
        rho2D[i] = (double*) malloc(NGy*sizeof(double));

        uSMonGird[i] = (Vec**) malloc(NGy*sizeof(Vec*));
  			newR[i] = (Vec**) malloc(NGy*sizeof(Vec*));
  			rhoOnGird[i] = (double**) malloc(NGy*sizeof(double*));

        for(int j = 0; j < NGy; j++){
  				uSMonGird[i][j] = (Vec*) malloc(NGz*sizeof(Vec));
  				newR[i][j] = (Vec*) malloc(NGz*sizeof(Vec));
  				rhoOnGird[i][j] = (double*) malloc(NGz*sizeof(double));
  			}
  		}



      cout << "loading data @ " << tCurrent+t_step << endl;
      loadData(tCurrent+t_step, pos2, posU2);

      for(int i = 0; i < NRPART; i++){
        uTot[i] = posU2[i] - posU[i];
      }

      //sprintf(fname, "uDisc-%d.dat",tCurrent);
      //StoreUField(t,fname,pos,uTot);

      cout << "making grid based displacement field" << endl;
	    makeGridBased(tCurrent, pos, uTot, newR, uSMonGird, rhoOnGird,0., NGx, NGy, NGz, delLoc);

      //sprintf(fname, "uCG-%d.dat",tCurrent);
      //storeGBVecForOvito(fname,t, newR, uSMonGird,rhoOnGird, NGx, NGy, NGz);
      cout << "projecting the fields onto a 2D grid" << endl;
      Make2DGrid(uSMonGird, v2D,rhoOnGird, rho2D, NGx, NGy, NGz);
      sprintf(fname, "uCG-2D-%d.dat",tCurrent);
      Store2DVec(t, fname,NGx, NGy, delLoc ,v2D, rho2D);
      cout << "calculating the shear strain on the grid" << endl;
      calShearStr2D(v2D, dUxy, dUyx, eps2D,NGx, NGy);
      cout << "storing the shear components" << endl;
      sprintf(fname, "strComponents-%d.dat",tCurrent);
      store2DShearComponents(t ,fname, dUxy, dUyx, eps2D, rho2D, NGx, NGy, delLoc );

      /*

      cout << "done with the griding " << endl;

      cout << "making a 2D projection" << endl;
      Make2DGrid(uSMonGird, v2D , NGx, NGy, NGz);
      //sprintf(fname, "Uz-%d.dat", tCurrent);

      CalGradUzx(v2D, dUzx, NGx, NGz);
      sprintf(fname, "DUzDx-%d.dat", tCurrent);
      Store2DScalar(tCurrent, fname,NGx,NGz,delLoc ,dUzx);
      CalGradUxz(v2D, dUxz,NGx,NGz);
      sprintf(fname, "DUxDz-%d.dat",tCurrent);
      Store2DScalar(tCurrent,fname,NGx,NGz,delLoc ,dUxz);
      for(int i = 0; i < NGx; i++){
        for(int k = 0; k < NGz; k++){
          eps2D[i][k]=(dUxz[i][k]+dUzx[i][k])*0.5;
        }
      }

      //sprintf(fname, "eps-%d.dat", tCurrent);
      //Store2DScalar(tCurrent,fname,NGx,NGz,delLoc,eps2D);

      IdentifySB(dUzx,actFlag, NGx,NGz);
      cout << "make 1D " << endl;
      MakeDUzx1D(actFlag,dUzx,dUzx1D,nr_dUzx1D,NGx,NGz);
      cout << "making dist. hist." << endl;
      calcHistDis(t, dUzx1D, nr_dUzx1D, histDis, NGx, NGz,delLoc,sizeHistDis);
      */

      for(int i = 0; i < NGx; i++){
        free(dUxy[i]);// = (double*)malloc(NGz*sizeof(double));
        free(dUyx[i]);// = (double*)malloc(NGz*sizeof(double));
        free(eps2D[i]);// = (double*)malloc(NGz*sizeof(double));
        free(actFlag[i]);// = (int*)malloc(NGz*sizeof(int));
        free(rho2D[i]);
      }
      free(dUyx);// = (double*)malloc(NGz*sizeof(double));
      free(dUxy);// = (double*)malloc(NGz*sizeof(double));
      free(eps2D);// = (double*)malloc(NGz*sizeof(double));
      free(actFlag);
      free(rho2D);


      free(dUzx1D);
      free(nr_dUzx1D);

      free(corrDu_sum);
      free(nrCorrDu_sum);

      free(corrDu);
      free(nrCorrDu);

      for(int i = 0; i < NGx; i++){
        for(int j = 0; j < NGy; j++){
          free(uSMonGird[i][j]);// = (Vec*) malloc(NGz*sizeof(Vec));
          free(newR[i][j]);// = (Vec*) malloc(NGz*sizeof(Vec));
          free(rhoOnGird[i][j]);
        }
        free(v2D[i]);// = (Vec*) malloc(NGz*sizeof(Vec));
        free(uNon[i]);
        free(uSMonGird[i]);// = (Vec*) malloc(NGz*sizeof(Vec));
        free(newR[i]);// = (Vec*) malloc(NGz*sizeof(Vec));
        free(rhoOnGird[i]);
      }

      free(v2D);// = (Vec*) malloc(NGz*sizeof(Vec));
      free(uNon);
      free(uSMonGird);// = (Vec*) malloc(NGz*sizeof(Vec));
      free(newR);// = (Vec*) malloc(NGz*sizeof(Vec));
      free(rhoOnGird);

    }
    sprintf(fname,"avgCorr.dat");
    //Save1DCorr(fname,corrDu_sum, nrCorrDu_sum,NGx, delLoc);
    free(pos);free(posU);//free(vel);

		free(pos2);free(posU2);//free(vel2);

    free(epsList);

    /*for(int i = 0; i < NRPART; i++){
      free(nblist[i]);
      free(Largenblist[i]);
    }

    free(nblist);
    free(Largenblist);*/


    end = clock();
    cout << "Time required for execution:"<< (double)(end-start)/CLOCKS_PER_SEC<< " seconds." << endl;

	  cout << "DONE." << endl;
    return 0;
  }
