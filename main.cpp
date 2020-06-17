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

string listInput[] = {"t_init", "t_step","num_times","NRPART","cut_off","dump_path","periodicity", "normalized_coordinate", "proj_plane","save_3D"};
const int NLIST_digit = 5;
const int NLIST = 10;
double incrStr;
int NRPART;
double LX;
double LY;
double LZ;

int t_first;
int num_times;
int t_step;
int t_delta;
string dump_path;
int dump_path_size;

int flag_t_first, flag_t_step, flag_num_times, flag_NRPART, flag_cutOff;
int flag_dump, flag_per[3], flag_normalization, flag_proj, flag_save_3D;


double intvl;
double gdot;

//const int MAX_NB=200;
const int MAX_NB_Large = 400;
double CUTOFF;
double BOX_CUTOFF;
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
	dr = p1 - p2;
  if(flag_per[0] == 1){
    dr.x -= ((int)(2.*dr.x*LXINV0)) * LX0 ;
  }
  if(flag_per[1] ==  1){
    dr.y -= ((int)(2.*dr.y*LYINV0)) * LY0 ;
  }
  if(flag_per[2] == 1){
    dr.z -= ((int)(2.*dr.z*LZINV0)) * LZ0 ;
  }

}
void apply_pbc(int time, Vec& mol)
{
  if(flag_per[0] == 1){
      mol.x = mol.x - ((int) floor((mol.x - XMIN) / LX)) * LX;
  }
  if(flag_per[1] == 1){
      mol.y = mol.y - ((int) floor((mol.y - YMIN) / LY)) * LY;
  }
  if(flag_per[2] == 1){
	   mol.z = mol.z - ((int) floor((mol.z - ZMIN) / LZ)) * LZ;
  }
}
void compare_assign(int caseNr, string val){
	//int flag_bondPath, flag_initData, flag_time_series, flag_dump;
	switch (caseNr) {
		case NLIST_digit:
    dump_path_size = val.length();
    dump_path = val;
		cout << "dump_path= " << dump_path<< endl;
		flag_dump = 1;
		break;
    case NLIST_digit+1:
    {
      int ind_flag = 0;
      for(int i = 0; i < val.length(); i++){
        if(val[i] != ' '){
          if(val[i] == 'p' || val[i] == 's' || val[i]=='f'){
            if(val[i] == 'p'){
              flag_per[ind_flag] = 1;
            }else if(val[i]=='s' || val[i]=='f'){
              flag_per[ind_flag]=2;
            }
            ind_flag++;
          }else{
            cout << "Error: "<< endl;
            cout << "the given input for the periodicity is not correct" << endl;
            cout << "the given values for each direction must be from {\'s\', \'p\',\'f\'}" << endl;
            cout << "given input: "<< val << endl;
            exit(16);
          }
        }
      }
      if(flag_per[0] == 1){
        cout << "periodic in x direction" << endl;
      }else{
        cout << "non-periodic in x direction" << endl;
      }
      if(flag_per[1] == 1){
        cout << "periodic in y direction" << endl;
      }else{
        cout << "non-periodic in y direction" << endl;
      }
      if(flag_per[2] == 1){
        cout << "periodic in z direction" << endl;
      }else{
        cout << "non-periodic in z direction" << endl;
      }
      break;
    }
		case NLIST_digit+2:
		if(val == "yes"){
			flag_normalization = 1;
		}else if(val == "YES"){
			flag_normalization = 1;
		}else if(val == "Yes"){
			flag_normalization = 1;
		}
		cout << "flag_normalization= " << flag_normalization << endl;
		break;
    case NLIST_digit+3:
    if(val == "xy" || val == "yx"){
      flag_proj = 1;
    }else if(val == "xz" || val == "zx"){
      flag_proj = 2;
    }else if(val == "zy"|| val == "yz"){
      flag_proj = 3;
    }else if(val == "none"){
      flag_proj = 4;
    }
    case NLIST_digit+4:
    if(val == "yes"){
			flag_save_3D =1;
		}else if(val == "YES"){
			flag_save_3D =1;
		}else if(val == "Yes"){
			flag_save_3D =1;
		}else if(val == "no"){
      flag_save_3D = 2;
    }
    //static int count= 0;
    //count++;
		cout << "flag_save_3D = " << flag_save_3D << endl;
		break;

	}
}
void compare_assign(int caseNr, double val){

	switch (caseNr) {
    case 0:
    t_first = (int) val;
    flag_t_first=1;
    cout << "t_first= " << t_first << endl;
    break;
    case 1:
    t_step = (int) val;
    flag_t_step = 1;
    cout << "t_step= " << t_step << endl;
    break;
		case 2:	//LX = val;
		num_times = (int)val;
    flag_num_times = 1;
		cout << "num_times= " << num_times << endl;
		break;
		case 3:
    NRPART = (int) val;
		cout << "NRPART= " << NRPART<< endl;
		flag_NRPART = 1;
		break;
    case 4:
    CUTOFF = val;
    BOX_CUTOFF=2.2*CUTOFF;
    cout << "cutoff: " << CUTOFF << endl;
    flag_cutOff = 1;
    break;
	}
}
void CheckInput(){
  int exit_flag= 0;
  if(flag_dump == 0){
    cout << "the dump file path is not set" << endl;
    exit_flag++;
  }
  if(flag_NRPART== 0){
    cout << "the number of particles is not set" << endl;
    exit_flag++;
  }
  if(flag_t_step== 0){
    cout << "the incremental timestep is not set" << endl;
    exit_flag++;
  }
  if(flag_t_first== 0){
    cout << "the initial timestep not set" << endl;
    exit_flag++;
  }
  if(flag_num_times== 0){
    cout << "the number of times is not set" << endl;
    exit_flag++;
  }
  if(flag_normalization == 0){
    cout << "the normalizitian flag is not set" << endl;
    exit_flag++;
  }
  if(flag_save_3D == 0){
    cout << "the preference for 3D save is not set" << endl;
    exit_flag++;
  }
  for(int i = 0; i < 3; i++){
    if(flag_per[i] == 0){
      cout << "periodicity is not set correctly" << endl;
      exit_flag++;
    }
  }
  if(flag_proj == 0){
    cout << "2D projection is not set correctly" << endl;
    exit_flag++;
  }
  if(flag_cutOff == 0){
    cout << "cutoff is not set correctly" << endl;
    exit_flag++;
  }
  if(exit_flag > 0){
    cout << "number of error with the input file : " << exit_flag << endl;
    exit(1);
  }
}
void readInput(char *fname){
  /*
  This function reads the input file and assign the needed parameters from the input file.
  */
	ifstream inFile(fname);
	if(!inFile){
		cout << "the initial input file could not be found" << endl;
		cout << "file name: "<< fname << endl;
		exit(1);
 	}

	string str;

	double val;
	int i = 0;//line number
	int counter = 0;// number of correct inputs
	size_t spac;
	string sub1;
	string sub2;

  flag_t_first = 0;flag_t_step=0; flag_num_times=0; flag_NRPART=0;
  flag_dump = 0;flag_normalization=0;flag_cutOff = 0;flag_save_3D = 0;
  flag_proj = 0;
  for(int i = 0; i < 3; i++)
    flag_per[i]=0;

	while(!inFile.eof()){
	 	getline(inFile, str);
		spac = str.find("=");
		sub1 = str.substr(0,spac);
		sub2 = str.substr(spac+1);
		int flag = 0;
		for(int j = 0; j < NLIST; j++){
			if(sub1.compare(listInput[j]) == 0){
				if(j < NLIST_digit){
					val = stod(sub2);
					compare_assign(j,val);
				}else{
					compare_assign(j,sub2);
				}
				flag = 1;
				counter++;
      }
		}
		if(flag == 0 && str.length() != 0){
			cout << "length: " << str.length() << endl;
			cout << "The input file at line " << i+1 << " is not defined well" << endl;
			cout << "Take a look!" << endl;
			cout << "str: " << str << endl;
			spac = str.find("=");
			sub1 = str.substr(0,spac);
			sub2 = str.substr(spac+1);
			cout << "sub1: " << sub1 << " sub2:" << sub2 << endl;
			cout << "spac: "<< spac << endl;
			exit(2);
		}
		i++;
	 }
   CheckInput();
}

void set_dumpFile_path(int t, char *arr){
  size_t spac;
  string sub1;
  string sub2;
  spac = dump_path.find("*");
  sub1 = dump_path.substr(0,spac);
  sub2 = dump_path.substr(spac+1);
  //cout << dump_path << endl;
  //cout << sub1 << endl;
  //cout << sub2 << endl;
  char arr2[100];
  char arr3[100];
  /*for(int i = 0; i < sub1.length(); i++){
    arr2[i] = sub1[i];
  }
  for(int i = 0; i < sub2.length(); i++){
    arr3[i] = sub2[i];
  }*/
  //cout << arr2 << endl;
  strcpy(arr2, sub1.c_str());
  strcpy(arr3, sub2.c_str());
  sprintf(arr, "%s%d%s",arr2,t,arr3);
  cout << arr << endl;
  int n = sub1.length()+sub2.length()+7;
  //cout << "size: " << n << endl;
  arr[n]='\0';
}

void loadData(int t, Vec* data, Vec* posU)
{
	int i, iniLine= 5, tempInd;
  string tempstr;
  //memcpy(tempstr,dump_path,dump_path_size);

  char fname[1200];
  set_dumpFile_path(t, fname);
	//sprintf(fname, "%s%d%s",sub1,sub2, t);
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
  double x, y, z, tmpD;
  inputFile.precision(15);
	for(i = 0; i < NRPART; i++){
		inputFile >> tempInd;
		inputFile >> typen[tempInd - 1] >>  x >> y >> z >> tmpD >> tmpD >> tmpD;
    if(i == 0) {
      cout << tempInd << " " << typen[tempInd - 1] << " " <<  x << " " << y << " " << z << endl;
    }

    if(flag_normalization == 1){
      data[tempInd - 1].x = x * LX + XMIN;
      data[tempInd - 1].y = y * LY + YMIN;
      data[tempInd - 1].z = z * LZ + ZMIN;
    }else{
      data[tempInd - 1].x = x ;
      data[tempInd - 1].y = y ;
      data[tempInd - 1].z = z ;
    }
		posU[tempInd - 1] = data[tempInd - 1];
		apply_pbc(t,data[tempInd - 1]);
	}
	inputFile.close();
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
              if(flag_per[0] == 1){
                if(ixn < 0){
                  ixn += NR_SUBX;
                }else if(ixn >= NR_SUBX){
                  ixn -= NR_SUBX;
                }
              }else{
                if(ixn < 0){
                  continue;
                }else if(ixn >= NR_SUBX){
                  continue;
                }
              }
              //iyn = iy_help + iy;
              if(flag_per[1] == 1){
                if(iyn < 0){
                  iyn += NR_SUBX;
                }else if(iyn >= NR_SUBY){
                  iyn -= NR_SUBX;
                }
              }else{
                if(iyn < 0){
                  continue;
                }else if(iyn >= NR_SUBY){
                  continue;
                }
              }
              if(flag_per[2] == 1){
                if(izn < 0){
                  izn += NR_SUBZ;
                  //continue;
                }else if(izn >= NR_SUBZ){
                  izn -= NR_SUBZ;
                  //continue;
                }
              }else{
                if(izn < 0){
                  continue;
                }else if(izn >= NR_SUBZ){
                  continue;
                }
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

void Store2DVec(int t, char* fname,int NG1, int NG2, double delLoc ,Vec** u, double **rho){
	ofstream oFile(fname);
	oFile.precision(15);
  if(flag_proj == 1){
    oFile << "ITEM: TIMESTEP\n"<< t<<"\nITEM: NUMBER OF ATOMS\n"<<NG1*NG2<<"\nITEM: BOX BOUNDS xy xz yz pp pp pp\n"<< XMIN<< " "<< XMAX <<" 0\n"<< YMIN<< " "<< YMAX<<" 0\n"<< -1 << " "<< 1  << " 0"<<endl;
  }else if(flag_proj == 2){
    oFile << "ITEM: TIMESTEP\n"<< t<<"\nITEM: NUMBER OF ATOMS\n"<<NG1*NG2<<"\nITEM: BOX BOUNDS xy xz yz pp pp pp\n"<< XMIN<< " "<< XMAX <<" 0\n"<< -1 << " "<< 1<<" 0\n"<< ZMIN<< " "<< ZMAX  << " 0"<<endl;
  }else if(flag_proj == 3){
    oFile << "ITEM: TIMESTEP\n"<< t<<"\nITEM: NUMBER OF ATOMS\n"<<NG1*NG2<<"\nITEM: BOX BOUNDS xy xz yz pp pp pp\n"<< -1 << " "<< 1 <<" 0\n"<< YMIN<< " "<< YMAX<<" 0\n"<< ZMIN<< " "<< ZMAX  << " 0"<<endl;
  }
	oFile << "ITEM: ATOMS id type x y z ux uy uz rho\n";
  int nr = 1;
  if(flag_proj== 1){
  	for(int i = 0; i < NG1; i++){
      for(int j = 0; j < NG2; j++){
          oFile << nr << " "<< 1 << " " << XMIN+(i+0.5)*delLoc << " "<< YMIN+(j+0.5)*delLoc<< " 0. " << u[i][j].x << " "<<  u[i][j].y << " " << u[i][j].z << " " << rho[i][j]<< endl;
          nr++;
      }
  	}
  }else if(flag_proj == 2){
    for(int i = 0; i < NG1; i++){
      for(int j = 0; j < NG2; j++){
          oFile << nr << " "<< 1 << " " << XMIN+(i+0.5)*delLoc << " 0. "<< ZMIN+(j+0.5)*delLoc<< " " << u[i][j].x << " "<<  u[i][j].y << " " << u[i][j].z << " " << rho[i][j]<< endl;
          nr++;
      }
  	}
  }else if(flag_proj == 3){
    for(int i = 0; i < NG1; i++){
      for(int j = 0; j < NG2; j++){
          oFile << nr << " "<< 1 << " 0. "  << YMIN+(i+0.5)*delLoc<< " " << ZMIN+(j+0.5)*delLoc<< " "<< u[i][j].x << " "<<  u[i][j].y << " " << u[i][j].z << " " << rho[i][j]<< endl;
          nr++;
      }
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
void storeGBVecCSVStyle(char *fname, int t, Vec ***r, Vec ***u, int Nx, int Ny, int Nz){
	ofstream oFile(fname);
	oFile.precision(15);
	oFile << "X,Y,Z,Ux,Uy,Uz\n";
	//int nr = 0;
	for(int i = 0; i < Nx; i++){
		for(int j = 0; j < Ny; j++){
			for(int k = 0; k < Nz; k++){
					oFile << r[i][j][k].x << ","<< r[i][j][k].y << ","<< r[i][j][k].z << ","<< u[i][j][k].x <<" "<< u[i][j][k].y <<" " << u[i][j][k].z << endl;//" " << UField[i].length()

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
	int nr = 1;
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
void storeGBTenForOvito(char *fname, int t, Vec ***r, p_ten2 ***eps, double ***rho, int Nx, int Ny, int Nz){
  ofstream oFile(fname);
	oFile.precision(15);
	oFile << "ITEM: TIMESTEP\n"<< t <<"\nITEM: NUMBER OF ATOMS\n"<<Nx*Ny*Nz<<"\nITEM: BOX BOUNDS xy xz yz pp pp pp\n"<< XMIN0<< " "<< XMAX0 <<" 0\n"<< YMIN0 << " "<< YMAX0 <<" 0\n"<< ZMIN0<< " "<< ZMAX0 <<" 0"<< endl;
	oFile << "ITEM: ATOMS id type x y z sxx sxy sxz syx syy syz szx szy szz rho\n";
	int nr = 1;
	for(int i = 0; i < Nx; i++){
		for(int j = 0; j < Ny; j++){
			for(int k = 0; k < Nz; k++){
					oFile << nr << " 1 " << r[i][j][k] <<  " " << eps[i][j][k]<< " "<<  rho[i][j][k] << endl;//" " << UField[i].length()
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
	int nr = 1;
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
void Make2DGrid(Vec ***u, Vec **v,double ***rho3D ,double **rho2D, int Nx, int Ny, int Nz){
  if(flag_proj == 1){
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
  }else if(flag_proj == 2){
    for(int i = 0; i < Nx; i++){
        for(int k = 0; k < Nz; k++){
          v[i][k].x = 0.;
          v[i][k].z = 0.;
          v[i][k].y = 0.;
          rho2D[i][k] = 0.;
          for(int j = 0; j < Ny; j++){
            v[i][k].x += u[i][j][k].x;
            v[i][k].z += u[i][j][k].z;
            v[i][k].y += u[i][j][k].y;
            rho2D[i][k] += rho3D[i][j][k];
          }
          v[i][k].x = v[i][k].x / (double)Ny;
          v[i][k].z = v[i][k].z / (double)Ny;
          v[i][k].y = v[i][k].y / (double)Ny;
          rho2D[i][k] = rho2D[i][k] / (double)Ny;
      }
    }
  }else if(flag_proj == 3){
    for(int j = 0; j < Ny; j++){
        for(int k = 0; k < Nz; k++){
          v[j][k].x = 0.;
          v[j][k].z = 0.;
          v[j][k].y = 0.;
          rho2D[j][k] = 0.;
          for(int i = 0; i < Nx; i++){
            v[j][k].x += u[i][j][k].x;
            v[j][k].z += u[i][j][k].z;
            v[j][k].y += u[i][j][k].y;
            rho2D[j][k] += rho3D[i][j][k];
          }
          v[j][k].x = v[j][k].x / (double)Nx;
          v[j][k].z = v[j][k].z / (double)Nx;
          v[j][k].y = v[j][k].y / (double)Nx;
          rho2D[j][k] = rho2D[j][k] / (double)Nx;
      }
    }
  }
}
void CalGradUyx(Vec **v, double **dUyx, int Nx, int Ny,double delLoc){
  //double delLoc = 1.;
  for(int i = 0; i < Nx; i++){
    int iMinus = i-1;
    int iPlus = i+1;
    if(flag_per[0]=1){
      if(iMinus == -1){
        iMinus += Nx;
      }else if(iPlus == Nx){
        iPlus -= Nx;
      }
      for(int j = 0; j < Ny; j++){
        dUyx[i][j] = 0.5/(delLoc)*(v[iPlus][j].y - v[iMinus][j].y);
      }
    }else{
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
}
void CalGradUxy(Vec **v, double **dUxy, int Nx, int Ny, double delLoc){
  //double delLoc = 1.;
  for(int j = 0; j < Ny; j++){
    int jMinus = j-1;
    int jPlus = j+1;
    if(flag_per[1] == 1){
      if(jMinus == -1){
        jMinus += Ny;
      }else if(jPlus == Ny){
        jPlus -= Ny;
      }
      for(int i = 0; i < Nx; i++){
        dUxy[i][j] = 0.5/(delLoc)*(v[i][jPlus].x - v[i][jMinus].x);
      }
    }else{
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
}
void CalGradUxz(Vec **v, double **dUxz, int Nx, int Nz, double delLoc){
  //double delLoc = 1.;
  for(int k = 0; k < Nz; k++){
    int kMinus = k-1;
    int kPlus = k+1;
    if(flag_per[2] == 1){
      if(kMinus == -1){
        kMinus += Nz;
      }else if(kPlus == Nz){
        kPlus -= Nz;
      }
      for(int i = 0; i < Nx; i++){
        dUxz[i][k] = 0.5/(delLoc)*(v[i][kPlus].x - v[i][kMinus].x);
      }
    }else{
      if(kMinus == -1){
        for(int i = 0; i < Nx; i++){
          dUxz[i][k] = 1./(delLoc)*(v[i][kPlus].x - v[i][k].x);
        }
      }else if(kPlus == Nz){
        for(int i = 0; i < Nx; i++){
          dUxz[i][k] = 1./(delLoc)*(v[i][k].x - v[i][kMinus].x);
        }
      }else{
        for(int i = 0; i < Nx; i++){
          dUxz[i][k] = 0.5/(delLoc)*(v[i][kPlus].x - v[i][kMinus].x);
        }
      }
    }
  }
}
void CalGradUzx(Vec **v, double **dUzx, int Nx, int Nz,double delLoc){
  //double delLoc = 1.;
  for(int i = 0; i < Nx; i++){
    int iMinus = i-1;
    int iPlus = i+1;
    if(flag_per[0]=1){
      if(iMinus == -1){
        iMinus += Nx;
      }else if(iPlus == Nx){
        iPlus -= Nx;
      }
      for(int k = 0; k < Nz; k++){
        dUzx[i][k] = 0.5/(delLoc)*(v[iPlus][k].z - v[iMinus][k].z);
      }
    }else{
      if(iMinus == -1){
        iMinus = 0;
        for(int k = 0; k < Nz; k++){
          dUzx[i][k] = 1./(delLoc)*(v[iPlus][k].z - v[iMinus][k].z);
        }
      }else if(iPlus == Nx){
        iPlus = Nx-1;
        for(int k = 0; k < Nz; k++){
          dUzx[i][k] = 1./(delLoc)*(v[iPlus][k].z - v[iMinus][k].z);
        }
      }else{
        for(int k = 0; k < Nz; k++){
          dUzx[i][k] = 0.5/(delLoc)*(v[iPlus][k].z - v[iMinus][k].z);
        }
      }
    }
  }
}
void CalGradUyz(Vec **v, double **dUyz, int Ny, int Nz, double delLoc){
  //double delLoc = 1.;
  for(int k = 0; k < Nz; k++){
    int kMinus = k-1;
    int kPlus = k+1;
    if(flag_per[2] == 1){
      if(kMinus == -1){
        kMinus += Nz;
      }else if(kPlus == Nz){
        kPlus -= Nz;
      }
      for(int j = 0; j < Ny; j++){
        dUyz[j][k] = 0.5/(delLoc)*(v[j][kPlus].y - v[j][kMinus].y);
      }
    }else{
      if(kMinus == -1){
        for(int j = 0; j < Ny; j++){
          dUyz[j][k] = 1./(delLoc)*(v[j][kPlus].y - v[j][k].y);
        }
      }else if(kPlus == Nz){
        for(int j = 0; j < Ny; j++){
          dUyz[j][k] = 1./(delLoc)*(v[j][k].y - v[j][kMinus].y);
        }
      }else{
        for(int j = 0; j < Ny; j++){
          dUyz[j][k] = 0.5/(delLoc)*(v[j][kPlus].y - v[j][kMinus].y);
        }
      }
    }
  }
}
void CalGradUzy(Vec **v, double **dUzy, int Ny, int Nz,double delLoc){
  //double delLoc = 1.;
  for(int j = 0; j < Ny; j++){
    int jMinus = j-1;
    int jPlus = j+1;
    if(flag_per[1]=1){
      if(jMinus == -1){
        jMinus += Ny;
      }else if(jPlus == Ny){
        jPlus -= Ny;
      }
      for(int k = 0; k < Nz; k++){
        dUzy[j][k] = 0.5/(delLoc)*(v[jPlus][k].z - v[jMinus][k].z);
      }
    }else{
      if(jMinus == -1){
        jMinus = 0;
        for(int k = 0; k < Nz; k++){
          dUzy[j][k] = 1./(delLoc)*(v[jPlus][k].z - v[jMinus][k].z);
        }
      }else if(jPlus == Ny){
        jPlus = Ny-1;
        for(int k = 0; k < Nz; k++){
          dUzy[j][k] = 1./(delLoc)*(v[jPlus][k].z - v[jMinus][k].z);
        }
      }else{
        for(int k = 0; k < Nz; k++){
          dUzy[j][k] = 0.5/(delLoc)*(v[jPlus][k].z - v[jMinus][k].z);
        }
      }
    }
  }
}
void calShearStr2D(Vec **v, double **dU12, double **dU21, double **eps,int N1, int N2, double delLoc){
  if(flag_proj == 1){
    cout << "calculating yx component" << endl;
    CalGradUyx(v,dU21, N1,N2, delLoc);
    cout << "calculating xy component" << endl;
    CalGradUxy(v,dU12, N1,N2, delLoc);
    cout << "summing up" << endl;
    for(int i = 0; i < N1; i++){
      for(int j = 0; j < N2; j++){
        eps[i][j] = 0.5 * (dU21[i][j]+ dU12[i][j]);
      }
    }
  }else if(flag_proj == 2){
    cout << "calculating zx component" << endl;
    CalGradUzx(v,dU21, N1,N2, delLoc);
    cout << "calculating xz component" << endl;
    CalGradUxz(v,dU12, N1,N2, delLoc);
    cout << "summing up" << endl;
    for(int i = 0; i < N1; i++){
      for(int k = 0; k < N2; k++){
        eps[i][k] = 0.5 * (dU21[i][k]+ dU12[i][k]);
      }
    }
  }else if(flag_proj == 3){
    cout << "calculating zy component" << endl;
    CalGradUzx(v,dU21, N1,N2, delLoc);
    cout << "calculating zy component" << endl;
    CalGradUxz(v,dU12, N1,N2, delLoc);
    cout << "summing up" << endl;
    for(int j = 0; j < N1; j++){
      for(int k = 0; k < N2; k++){
        eps[j][k] = 0.5 * (dU21[j][k]+ dU12[j][k]);
      }
    }
  }
}
void cal3DGradX(Vec ***uSM, p_ten2 ***DU, int Nx, int Ny, int Nz, double delLoc){
  Vec du;
  for(int i = 0; i < Nx; i++){
    int iMinus = i - 1;
    int iPlus = i + 1;
    if(flag_per[0] == 1){
      if(iMinus == -1){
        iMinus += Nx;
      }else if(iPlus == Nx){
        iPlus -= Nx;
      }
      for(int j = 0; j < Ny; j++){
        for(int k = 0; k < Nz; k++){
          du = uSM[iPlus][j][k]-uSM[iMinus][j][k];
          DU[i][j][k].xx = 0.5*du.x/delLoc;
          DU[i][j][k].yx = 0.5*du.y/delLoc;
          DU[i][j][k].zx = 0.5*du.z/delLoc;
        }
      }
    }else{
      if(iMinus == -1){
        iMinus = i;
        for(int j = 0; j < Ny; j++){
          for(int k = 0; k < Nz; k++){
            du = uSM[iPlus][j][k]-uSM[iMinus][j][k];
            DU[i][j][k].xx = du.x/delLoc;
            DU[i][j][k].yx = du.y/delLoc;
            DU[i][j][k].zx = du.z/delLoc;
          }
        }
      }else if(iPlus == Nx){
        iPlus = i;
        for(int j = 0; j < Ny; j++){
          for(int k = 0; k < Nz; k++){
            du = uSM[iPlus][j][k]-uSM[iMinus][j][k];
            DU[i][j][k].xx = du.x/delLoc;
            DU[i][j][k].yx = du.y/delLoc;
            DU[i][j][k].zx = du.z/delLoc;
          }
        }
      }else{
        for(int j = 0; j < Ny; j++){
          for(int k = 0; k < Nz; k++){
            du = uSM[iPlus][j][k]-uSM[iMinus][j][k];
            DU[i][j][k].xx = 0.5*du.x/delLoc;
            DU[i][j][k].yx = 0.5*du.y/delLoc;
            DU[i][j][k].zx = 0.5*du.z/delLoc;
          }
        }
      }
    }
  }
}
void cal3DGradY(Vec ***uSM, p_ten2 ***DU, int Nx, int Ny, int Nz, double delLoc){
  Vec du;
  for(int j = 0; j < Ny; j++){
    int jMinus = j - 1;
    int jPlus = j + 1;
    if(flag_per[1] == 1){
      if(jMinus == -1){
        jMinus += Ny;
      }else if(jPlus == Ny){
        jPlus -= Ny;
      }
      for(int i = 0; i < Nx; i++){
        for(int k = 0; k < Nz; k++){
          du = uSM[i][jPlus][k]-uSM[i][jMinus][k];
          DU[i][j][k].xy = 0.5*du.x/delLoc;
          DU[i][j][k].yy = 0.5*du.y/delLoc;
          DU[i][j][k].zy = 0.5*du.z/delLoc;
        }
      }
    }else{
      if(jMinus == -1){
        jMinus = j;
        for(int i = 0; i < Nx; i++){
          for(int k = 0; k < Nz; k++){
            du = uSM[i][jPlus][k]-uSM[i][jMinus][k];
            DU[i][j][k].xy = du.x/delLoc;
            DU[i][j][k].yy = du.y/delLoc;
            DU[i][j][k].zy = du.z/delLoc;
          }
        }
      }else if(jPlus == Ny){
        jPlus = j;
        for(int i = 0; i < Nx; i++){
          for(int k = 0; k < Nz; k++){
            du = uSM[i][jPlus][k]-uSM[i][jMinus][k];
            DU[i][j][k].xy = du.x/delLoc;
            DU[i][j][k].yy = du.y/delLoc;
            DU[i][j][k].zy = du.z/delLoc;
          }
        }
      }else{
        for(int i = 0; i < Nx; i++){
          for(int k = 0; k < Nz; k++){
            du = uSM[i][jPlus][k]-uSM[i][jMinus][k];
            DU[i][j][k].xy = 0.5*du.x/delLoc;
            DU[i][j][k].yy = 0.5*du.y/delLoc;
            DU[i][j][k].zy = 0.5*du.z/delLoc;
          }
        }
      }
    }
  }
}
void cal3DGradZ(Vec ***uSM, p_ten2 ***DU, int Nx, int Ny, int Nz, double delLoc){
  Vec du;
  for(int k = 0; k < Nz; k++){
    int kMinus = k - 1;
    int kPlus = k + 1;
    if(flag_per[2] == 1){
      if(kMinus == -1){
        kMinus += Nz;
      }else if(kPlus == Nz){
        kPlus -= Nz;
      }
      for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
          du = uSM[i][j][kPlus]-uSM[i][j][kMinus];
          DU[i][j][k].xz = 0.5*du.x/delLoc;
          DU[i][j][k].yz = 0.5*du.y/delLoc;
          DU[i][j][k].zz = 0.5*du.z/delLoc;
        }
      }
    }else{
      if(kMinus == -1){
        kMinus = k;
        for(int i = 0; i < Nx; i++){
          for(int j = 0; k < Ny; j++){
            du = uSM[i][j][kPlus]-uSM[i][j][kMinus];
            DU[i][j][k].xz = du.x/delLoc;
            DU[i][j][k].yz = du.y/delLoc;
            DU[i][j][k].zz = du.z/delLoc;
          }
        }
      }else if(kPlus == Nz){
        kPlus = k;
        for(int i = 0; i < Nx; i++){
          for(int j = 0; j < Ny; j++){
            du = uSM[i][j][kPlus]-uSM[i][j][kMinus];
            DU[i][j][k].xz = du.x/delLoc;
            DU[i][j][k].yz = du.y/delLoc;
            DU[i][j][k].zz = du.z/delLoc;
          }
        }
      }else{
        for(int i = 0; i < Nx; i++){
          for(int j = 0; j < Ny; j++){
            du = uSM[i][j][kPlus]-uSM[i][j][kMinus];
            DU[i][j][k].xz = 0.5*du.x/delLoc;
            DU[i][j][k].yz = 0.5*du.y/delLoc;
            DU[i][j][k].zz = 0.5*du.z/delLoc;
          }
        }
      }
    }
  }
}
void cal3DStrain(Vec ***uSM, p_ten2 ***eps, int Nx, int Ny, int Nz, double delLoc){
    p_ten2 ***DU;
    DU = (p_ten2***) malloc(Nx*sizeof(p_ten2**));
    for(int i = 0; i < Nx; i++){
      DU[i] = (p_ten2**) malloc(Ny*sizeof(p_ten2*));
      for(int j = 0; j < Ny; j++){
        DU[i][j] = (p_ten2*) malloc(Nz*sizeof(p_ten2));
      }
    }
    cal3DGradX(uSM, DU, Nx, Ny, Nz, delLoc);
    cal3DGradY(uSM, DU, Nx, Ny, Nz, delLoc);
    cal3DGradZ(uSM, DU, Nx, Ny, Nz, delLoc);
    for(int i = 0; i < Nx; i++){
      for(int j = 0; j < Ny; j++){
        for(int k = 0; k < Nz; k++){
          eps[i][j][k]= (DU[i][j][k]+transpos(DU[i][j][k]))*0.5;
        }
      }
    }
    for(int i = 0; i < Nx; i++){
      for(int j = 0; j < Ny; j++){
        free(DU[i][j]);
      }
      free(DU[i]);
    }
    free(DU);
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
int main(int argc, char **argv)
{
	clock_t start, end;
	start = clock();
  if(argc !=2){
    	cout << "the path to the input file is missing from the given argument!!!" << endl;
      return 1;

  }
  char input_path[1000];
  strcpy(input_path, argv[1]);
  cout << "the input path: " << input_path << endl;
  readInput(input_path);
  cout << "NRPART: "<< NRPART << endl;
  char fname[1000];
  //load_dump_path(fname, "paths.txt");

  typen = (int*) malloc(NRPART*sizeof(int));

  //char fname[100];

  Vec *pos,*posU, *pos2,*posU2,*uTot;
  int **boxInd, **boxInd2;
  double *listD2;
  p_ten2 ***eps;

  uTot = (Vec*) malloc(NRPART*sizeof(Vec));
  pos = (Vec*) malloc(NRPART*sizeof(Vec));
  posU = (Vec*) malloc(NRPART*sizeof(Vec));

  pos2 = (Vec*) malloc(NRPART*sizeof(Vec));
  posU2 = (Vec*) malloc(NRPART*sizeof(Vec));

	int NGx,NGy,NGz;
  double delLoc = 3./4.*CUTOFF;

	Vec ***uSMonGird, ***newR;
	double ***rhoOnGird, **rho2D;
  Vec **v2D;
  double **dU21,**dU12, **eps2D;

  //Vec **uNon;
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
    if(flag_proj != 4){
      if(flag_proj == 1){
        dU12 = (double**)malloc(NGx*sizeof(double*));
        dU21 = (double**)malloc(NGx*sizeof(double*));
        eps2D = (double**)malloc(NGx*sizeof(double*));
        v2D = (Vec**) malloc(NGx*sizeof(Vec*));
        rho2D = (double**) malloc(NGx*sizeof(double*));
        for(int i = 0; i < NGx; i++){
          v2D[i] = (Vec*) malloc(NGy*sizeof(Vec));
          rho2D[i] = (double*) malloc(NGy*sizeof(double));
          dU12[i] = (double*)malloc(NGy*sizeof(double));
          dU21[i] = (double*)malloc(NGy*sizeof(double));
          eps2D[i] = (double*)malloc(NGy*sizeof(double));
        }
      }else if(flag_proj == 2){
        dU12 = (double**)malloc(NGx*sizeof(double*));
        dU21 = (double**)malloc(NGx*sizeof(double*));
        eps2D = (double**)malloc(NGx*sizeof(double*));
        v2D = (Vec**) malloc(NGx*sizeof(Vec*));
        rho2D = (double**) malloc(NGx*sizeof(double*));

        for(int i = 0; i < NGx; i++){
          dU12[i] = (double*)malloc(NGz*sizeof(double));
          dU21[i] = (double*)malloc(NGz*sizeof(double));
          eps2D[i] = (double*)malloc(NGz*sizeof(double));
          v2D[i] = (Vec*) malloc(NGz*sizeof(Vec));
          rho2D[i] = (double*) malloc(NGz*sizeof(double));
        }
      }else if(flag_proj == 3){
        dU12 = (double**)malloc(NGy*sizeof(double*));
        dU21 = (double**)malloc(NGy*sizeof(double*));
        eps2D = (double**)malloc(NGy*sizeof(double*));
        v2D = (Vec**) malloc(NGy*sizeof(Vec*));
        rho2D = (double**) malloc(NGy*sizeof(double*));

        for(int j = 0; j < NGy; j++){
          dU12[j] = (double*)malloc(NGz*sizeof(double));
          dU21[j] = (double*)malloc(NGz*sizeof(double));
          eps2D[j] = (double*)malloc(NGz*sizeof(double));
          v2D[j] = (Vec*) malloc(NGz*sizeof(Vec));
          rho2D[j] = (double*) malloc(NGz*sizeof(double));
        }
      }
    }

    uSMonGird= (Vec***) malloc(NGx*sizeof(Vec**));
		newR = (Vec***) malloc(NGx*sizeof(Vec**));
		rhoOnGird = (double***) malloc(NGx*sizeof(double**));
    eps = (p_ten2***) malloc(NGx*sizeof(p_ten2**));

		for(int i = 0; i < NGx; i++){
      uSMonGird[i] = (Vec**) malloc(NGy*sizeof(Vec*));
			newR[i] = (Vec**) malloc(NGy*sizeof(Vec*));
			rhoOnGird[i] = (double**) malloc(NGy*sizeof(double*));
      eps[i] = (p_ten2**) malloc(NGy*sizeof(p_ten2*));

      for(int j = 0; j < NGy; j++){
				uSMonGird[i][j] = (Vec*) malloc(NGz*sizeof(Vec));
				newR[i][j] = (Vec*) malloc(NGz*sizeof(Vec));
				rhoOnGird[i][j] = (double*) malloc(NGz*sizeof(double));
        eps[i][j] = (p_ten2*) malloc(NGz*sizeof(p_ten2));
			}
		}

    cout << "loading data @ " << tCurrent+t_step << endl;

    loadData(tCurrent+t_step, pos2, posU2);

    for(int i = 0; i < NRPART; i++){
      uTot[i] = posU2[i] - posU[i];
    }

    cout << "making grid based displacement field" << endl;
    makeGridBased(tCurrent, pos, uTot, newR, uSMonGird, rhoOnGird,0., NGx, NGy, NGz, delLoc);

    cout << "calculating strain on the 3D grid" << endl;
    cal3DStrain(uSMonGird, eps, NGx, NGy,NGz, delLoc);

    if(flag_save_3D == 1){
      cout << "storing 3D displacement" << endl;
      sprintf(fname, "3D_displ_onGrid_%d.dat", tCurrent);
      storeGBVecForOvito(fname, t, newR, uSMonGird, rhoOnGird, NGx, NGy,NGz);
      cout << "saving the 3D strain" << endl;
      sprintf(fname, "3D_strain_onGrid_%d.dat", tCurrent);
      storeGBTenForOvito(fname,t, newR, eps, rhoOnGird,NGx, NGy,NGz);
    }

    if(flag_proj != 4){
      if(flag_proj == 1){
        cout << "projecting the fields onto a 2D grid" << endl;
        Make2DGrid(uSMonGird, v2D,rhoOnGird, rho2D, NGx, NGy, NGz);
        sprintf(fname, "uCG-2D-%d.dat",tCurrent);
        Store2DVec(t, fname,NGx, NGy, delLoc ,v2D, rho2D);
        cout << "calculating the shear strain on the grid" << endl;
        calShearStr2D(v2D, dU12, dU21, eps2D,NGx, NGy, delLoc);
        cout << "storing the shear components" << endl;
        sprintf(fname, "strComponents-%d.dat",tCurrent);
        store2DShearComponents(t ,fname, dU12, dU21, eps2D, rho2D, NGx, NGy, delLoc );
      }else if(flag_proj == 2){
        cout << "projecting the fields onto a 2D grid" << endl;
        Make2DGrid(uSMonGird, v2D,rhoOnGird, rho2D, NGx, NGy, NGz);
        sprintf(fname, "uCG-2D-%d.dat",tCurrent);
        Store2DVec(t, fname,NGx, NGz, delLoc ,v2D, rho2D);
        cout << "calculating the shear strain on the grid" << endl;
        calShearStr2D(v2D, dU12, dU21, eps2D,NGx, NGz, delLoc);
        cout << "storing the shear components" << endl;
        sprintf(fname, "strComponents-%d.dat",tCurrent);
        store2DShearComponents(t ,fname, dU12, dU21, eps2D, rho2D, NGx, NGz, delLoc );
      }else if(flag_proj == 3){
        cout << "projecting the fields onto a 2D grid" << endl;
        Make2DGrid(uSMonGird, v2D,rhoOnGird, rho2D, NGx, NGy, NGz);
        sprintf(fname, "uCG-2D-%d.dat",tCurrent);
        Store2DVec(t, fname,NGy, NGz, delLoc ,v2D, rho2D);
        cout << "calculating the shear strain on the grid" << endl;
        calShearStr2D(v2D, dU12, dU21, eps2D,NGy, NGz, delLoc);
        cout << "storing the shear components" << endl;
        sprintf(fname, "strComponents-%d.dat",tCurrent);
        store2DShearComponents(t ,fname, dU12, dU21, eps2D, rho2D, NGy, NGz, delLoc );
      }
    }

    if(flag_proj != 4){
      if(flag_proj == 1 || flag_proj == 2){
          for(int i = 0; i < NGx; i++){
            free(dU12[i]);
            free(dU21[i]);
            free(eps2D[i]);
            free(rho2D[i]);
            free(v2D[i]);
          }
          free(dU21);
          free(dU12);
          free(eps2D);
          free(rho2D);
          free(v2D);
      }else if(flag_proj == 3){
        for(int j = 0; j < NGy; j++){
          free(dU12[j]);
          free(dU21[j]);
          free(eps2D[j]);
          free(rho2D[j]);
          free(v2D[j]);
        }
        free(dU21);
        free(dU12);
        free(eps2D);
        free(rho2D);
        free(v2D);
      }
    }

    for(int i = 0; i < NGx; i++){
      for(int j = 0; j < NGy; j++){
        free(uSMonGird[i][j]);
        free(newR[i][j]);
        free(rhoOnGird[i][j]);
        free(eps[i][j]);
      }

      free(uSMonGird[i]);
      free(newR[i]);
      free(rhoOnGird[i]);
      free(eps[i]);
    }

    free(uSMonGird);
    free(newR);
    free(rhoOnGird);
    free(eps);
  }
  free(pos);free(posU);
	free(pos2);free(posU2);
  free(uTot);

  end = clock();
  cout << "Time required for execution:"<< (double)(end-start)/CLOCKS_PER_SEC<< " seconds." << endl;

  cout << "DONE." << endl;
  return 0;
}
