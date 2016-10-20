//copyright by ArthurZhou @UMich&Fudan&HUST 
//1989-2006 
#include <stdlib.h>
#include <iostream>
#include <string>
#include <string.h>
#include <fstream>
#include <vector>
#include <cmath>

#include <sstream>
#include <algorithm>
#include <functional>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

//median
double mid(double *array, int N){
	int num;
	double middle;
	if(N%2==0){
		double num_value1, num_value2;
		num=N/2;
		middle=0;
		for(int i=0;i!=N;i++){
			int rank=1;
			int flag=0;
			for(int j=0;j!=N;j++){
				if(array[i]>array[j]){
					rank++;
				}
				else if(array[i]==array[j]&i!=j){
					flag++;
				}
			}
			if(rank==num){
				num_value1=array[i];
			}
			else if(rank<num&&(rank+flag)>=num){
				num_value1=array[i];
			}
			if(rank==num+1){
				num_value2=array[i];
			}
			else if(rank<num+1&&(rank+flag)>=num+1){
				num_value2=array[i];
			}
		}
		middle=(num_value1+num_value2)/2;
	}
	else if(N%2==1){
		num=(N+1)/2;
		middle=0;
		for(int i=0;i!=N;i++){
			int rank=1;
			int flag=0;
			
			//cout<<array[i]<<endl;
			
			for(int j=0;j!=N;j++){
				if(array[i]>array[j]){
					rank++;
				}
				else if(array[i]==array[j]&i!=j){
					flag++;
				}
			}
			
			if(rank==num){
				middle=array[i];
			}
			else if(rank<num&&(rank+flag)>=num){
				middle=array[i];
			}
			
		}
		
		
	}
	
	return middle;
	
	
} 






int main(int argc, char *argv[]){


//parameters and files    
	int bin;
	bin = 1000;
	int a,b;
	a=b=0;
	
	ifstream file2;
	ifstream file1;
	ofstream file3;
	//ofstream file4;
	
	file3.open("SFC_output.txt");
	//file4.open("SFC_pvalue.txt");
	
	for(int i=1;i=!argc;i++){
		if(strncmp(argv[i],"-b",2)==0){
			bin=atoi(argv[i+1]);
		}
		if(strncmp(argv[i],"-case",2)==0){
			a=atoi(argv[i+1]);
		}
		if(strncmp(argv[i],"-control",2)==0){
			b=atoi(argv[i+1]);
		}
		if(strncmp(argv[i],"-i",2)==0){
			file1.open(argv[i+1],ios::in | ios::binary);
			if(!file1.is_open()){
				cout<<"CANNOT OPEN INPUT FILE."<<endl;
				continue;
			}
			file2.open(argv[i+1],ios::in | ios::binary);
		}
		
	} 
	
	if(a==0||b==0){
		cout<<"Number of case/control can not be zero."<<endl;
		//break;
	}	
	
	int per_shift;
	int N_out;
	double k;
	int out;
	double shift;
	
//parameter define	
    string input;
    int line;
    for(int i=0;!file2.eof();i++){
    	getline(file2,input);
    	line=i; 
    }
    
    cout<<"There are "<<line<<" samples."<<endl;
    
    
    double **sample1;
	sample1=new double*[line];
	for(int i=0;i!=line;i++)sample1[i]=new double[a];
	
	double **sample2;
	sample2=new double*[line];
	for(int i=0;i!=line;i++)sample2[i]=new double[b];
	
    string **info;
	info=new string*[line];
	for(int i=0;i!=line;i++)info[i]=new string[2];
   
    file2.close();
    file2.clear();
   
    for(int i=0;i!=line;i++){
  		for(int j=0;j!=a;j++){
			file1>>sample1[i][j];
		}
		for(int j=0;j!=b;j++){
			file1>>sample2[i][j];
		}
    }
    
    
//SFC    

//ranking	
	double **cal;
	cal=new double*[line];
	for(int i=0;i!=line;i++) cal[i]=new double[6];
    
    int *rank;
	rank=new int[line];
    double *sample_average;
	sample_average=new double[line];
	
	for(int i=0;i!=line;i++){
		sample_average[i]=0;
		rank[i]=1; 
	}
	
	
	for(int i=0;i!=line;i++){
		for(int j=0;j!=a;j++){
			
			sample_average[i]=sample1[i][j]+sample_average[i];
		}
		for(int j=0;j!=b;j++){
			
			sample_average[i]=sample2[i][j]+sample_average[i];
		}
		
		sample_average[i]=sample_average[i]/(double)(a+b); 
	}	
    
    
    for(int i=0;i!=line;i++){
  		for(int j=0;j!=line;j++){
  			if(sample_average[j]<=sample_average[i]&&i!=j){
  				rank[i]++;
  			}
  		}
  	}
    
    for(int i=0;i!=line;i++){
  			for(int k=0;k!=line;k++){
  				if(rank[i]==rank[k]&&i!=k){
  					rank[k]--;
  				}
  			
  		}
  	} 
    
    
//median    
    for(int i=0;i!=line;i++){		
		cal[i][0]= mid(sample1[i],a)-mid(sample2[i],b);   	
    	cal[i][1]=cal[i][0]*cal[i][0]; 
    }
  
//SFC  
	
	int flank;
	flank=bin/2;
	
    for(int i=1;i!=(line+1);i++){	
    	for(int j=0;j!=line;j++){
    		if(rank[j]==i&&i<(flank+1)){
    			double *cal_mid;
				cal_mid=new double[i+flank];
    			for(int k=1;k!=(i+flank+1);k++){
    				int loc;
    				for(int l=0;l!=line;l++){
    					if(rank[l]==k){
    						loc=l;
    					}
    				}
    				cal_mid[k-1]=cal[loc][1];
    			}
    			cal[j][2]=mid(cal_mid,i+flank);
    			cal[j][2]=cal[j][2]/0.455;
    			cal[j][3]=sqrt(cal[j][2]);
				cal[j][4]=cal[j][0]/cal[j][3]; 
    			delete []cal_mid;
  			}
    		else if(rank[j]==i&&(line-i)<(flank+1)){
    			double *cal_mid;
				cal_mid=new double[line-i+flank+1];
    			for(int k=i-flank;k!=(line+1);k++){
    				int loc;
    				for(int l=0;l!=line;l++){
    					if(rank[l]==k){
    						loc=l;
    					}
    				}
    				cal_mid[k-i+flank]=cal[loc][1];
    			}
				cal[j][2]=mid(cal_mid,line-i+flank+1);
    			cal[j][2]=cal[j][2]/0.455;
    			cal[j][3]=sqrt(cal[j][2]);
				cal[j][4]=cal[j][0]/cal[j][3]; 
    			delete []cal_mid;
    		}
    		else if(rank[j]==i){
    			double *cal_mid;
				cal_mid=new double[flank+1];
    			for(int k=i-flank;k!=i+flank+1;k++){
    				int loc;
    				for(int l=0;l!=line;l++){
    					if(rank[l]==k){
    						loc=l;
    					}
    				}
    				cal_mid[k-i+flank]=cal[loc][1];
    			}
    			cal[j][2]=mid(cal_mid,flank+1);
    			cal[j][2]=cal[j][2]/0.455;
    			cal[j][3]=sqrt(cal[j][2]);
				cal[j][4]=cal[j][0]/cal[j][3]; 
    			delete []cal_mid;
    		}
    	}
    }
    
	for(int i=0;i!=line;i++){
		file3<<cal[i][4]<<endl;
  	}
    
    
    
    file1.close();
    file2.close();
    file3.close();
    //getchar(); 
}
