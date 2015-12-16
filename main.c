#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX_QUEUE 10000

double Px[128*128],Py[128*128],Pz[128*128];
unsigned char RangeImage[128*128];
unsigned char CopyImage[128*128];
int seed_matrix[128*128];//matrix stores label flags for seeding
double Sx[128*128],Sy[128*128],Sz[128*128];//surface normal matrices
int count;
double averageX,averageY,averageZ;
int threshold=150;//intensity for thresholding
int distance=2;//distance between pixels for cross product
double angle_threshold=0.9;//angle threshold in radians for joining pixel to region
int nextLevelIntensity=0;

int main()
{
void convertTo3D();
void RegionGrow(int,int);
double getAngle(double,double,double,double,double,double);
FILE *fpt;
char header[320];
char filename[30];
int    ROWS,COLS,BYTES;
int    r,c,index,i,j,k=1,flag;
double ax,bx,cx,ay,by,cy,az,bz,cz,vx1,vx2,vy1,vy2,vz1,vz2;

printf("Enter 128x128 input range image name: ");
gets(filename);

printf("\nReading input image...\n");
if ((fpt=fopen(filename,"rb")) == NULL)
{
  printf("Unable to open %s for reading\n",filename);
  exit(0);
}
fscanf(fpt,"%s %d %d %d\n",header,&COLS,&ROWS,&BYTES);
if (strcmp(header,"P5") != 0  ||  BYTES != 255 || ROWS!=128 || COLS!=128)
{
  printf("Not a greyscale 8-bit 128x128 PPM image\n");
  exit(0);
}

header[0]=fgetc(fpt);
fread(RangeImage,1,128*128,fpt);
fclose(fpt);
printf("Done.\n\n");

printf("Thresholding to remove background..\n");
for(r=0;r<128;r++)
    for(c=0;c<128;c++)
    {
        index=r*128+c;
        //Intensity of 150 as threshold works well, but also used some location based manual
        //fixes to remove noise i.e r<40,(r<60&&c<50),(r<60&&c>100) regions are taken as background
        if(r<40 || (r<60&&c<50) || (r<60&&c>100) || RangeImage[index]>threshold)
            RangeImage[index]=255;
    }
printf("Done..\n\n");

for(i=0;i<128*128;i++)
    CopyImage[i]=RangeImage[i];

fpt=fopen("thresh.ppm","w");
fprintf(fpt,"P5 %d %d 255\n",128,128);
fwrite(RangeImage,128*128,1,fpt);
fclose(fpt);
printf("Saved thresholded image as thresh.ppm\n\n");

printf("Converting range image to 3D coordinates..");
convertTo3D();
printf("\nDone..\n\n");

fpt=fopen("x.txt","w");
for(r=0;r<128;r++)
{
    for(c=0;c<128;c++)
        fprintf(fpt,"%3.0f    ",Px[r*128+c]);
    fprintf(fpt,"\n");
}
fclose(fpt);
fpt=fopen("y.txt","w");
for(r=0;r<128;r++)
{
    for(c=0;c<128;c++)
        fprintf(fpt,"%3.0f    ",Py[r*128+c]);
    fprintf(fpt,"\n");
}
fclose(fpt);
fpt=fopen("z.txt","w");
for(r=0;r<128;r++)
{
    for(c=0;c<128;c++)
        fprintf(fpt,"%3.0f    ",Pz[r*128+c]);
    fprintf(fpt,"\n");
}
fclose(fpt);
printf("Saved 3d coordinates as x.txt,y.txt,z.txt\n");


//initialize labels for seeding matrix
for(r=0;r<128;r++)
{
    for(c=0;c<128;c++)
    {
        index=r*128+c;
        if(RangeImage[index]==255)
            seed_matrix[index]=1;//label
        else
            seed_matrix[index]=0;//no label
    }
}


//calculate surface normals for needle map
for(r=0;r<128-distance;r++)
{
    for(c=0;c<128-distance;c++)
    {
        ax=Px[r*128+c];//current point 'a'
        bx=Px[r*128+c+distance];//point 'b' at distance 5 to the right of 'a'
        cx=Px[(r+distance)*128+c];//point 'c' at distance 5 below 'a'
        ay=Py[r*128+c];
        by=Py[r*128+c+distance];
        cy=Py[(r+distance)*128+c];
        az=Pz[r*128+c];
        bz=Pz[r*128+c+distance];
        cz=Pz[(r+distance)*128+c];
        vx1=bx-ax;//surface vector ab(x component)
        vx2=cx-ax;//surface vector ac(x component)
        vy1=by-ay;
        vy2=cy-ay;
        vz1=bz-az;
        vz2=cz-az;
        Sx[r*128+c]=vy1*vz2-vz1*vy2;//cross products of surface vectors to get surface normal
        Sy[r*128+c]=vz1*vx2-vx1*vz2;
        Sz[r*128+c]=vx1*vy2-vy1*vx2;
    }
}

fpt=fopen("surfaceX.txt","w");
for(r=0;r<128;r++)
{
    for(c=0;c<128;c++)
        fprintf(fpt,"%5.0f    ",Sx[r*128+c]);
    fprintf(fpt,"\n");
}
fclose(fpt);
fpt=fopen("surfaceY.txt","w");
for(r=0;r<128;r++)
{
    for(c=0;c<128;c++)
        fprintf(fpt,"%5.0f    ",Sy[r*128+c]);
    fprintf(fpt,"\n");
}
fclose(fpt);
fpt=fopen("surfaceZ.txt","w");
for(r=0;r<128;r++)
{
    for(c=0;c<128;c++)
        fprintf(fpt,"%5.0f    ",Sz[r*128+c]);
    fprintf(fpt,"\n");
}
fclose(fpt);
printf("Saved surface normals as surfaceX.txt,surfaceY.txt,surfaceZ.txt\n");

for(r=2;r<128-2-distance;r++)
    for(c=2;c<128-2-distance;c++)
    {
        flag=0;
        for(i=-2;i<2;i++)
            for(j=-2;j<2;j++)
                if(seed_matrix[(r+i)*128+(c+j)]==1)
                    flag=1;
        if(flag==0)//then this is seed pixel
        {
            RegionGrow(r,c);
            printf("\nTotal pixels in region %d are %d",k++,count);
            printf("\nAverage orientation vector is x:%4.2f, y:%4.2f, z:%4.2f\n",averageX,averageY,averageZ);
        }
    }
    printf("\n\nTotal regions are %d + 1 background wall",k-1);

    fpt=fopen("segment.ppm","w");
    fprintf(fpt,"P5 %d %d 255\n",128,128);
    fwrite(CopyImage,128*128,1,fpt);
    fclose(fpt);
    printf("\n\nSegmented image saved as segment.ppm\n");
    printf("Summary of parameters -\nIntensity threshold - 150\n");
    printf("Distance b/w pixels for cross product - 2\nAngle Threshold - 0.9 radians\n\n");

return 0;
}

double getAngle(double x1,double y1,double z1,double x2,double y2,double z2)
{
    double dotProduct=x1*x2+y1*y2+z1*z2;
    double magnitude1=sqrt(x1*x1+y1*y1+z1*z1);
    double magnitude2=sqrt(x2*x2+y2*y2+z2*z2);
    double theta=acos(dotProduct/(magnitude1*magnitude2));
    return theta;
}

void RegionGrow(int r,int c)
{
int	queue[MAX_QUEUE],qh,qt;
int newR,newC,r2,c2;
double diffAngle;
count=1;
seed_matrix[r*128+c]=1;

//initial average
averageX=Sx[r*128+c];
averageY=Sy[r*128+c];
averageZ=Sz[r*128+c];

queue[0]=r*128+c;
qh=1;	/* queue head */
qt=0;	/* queue tail */
while (qt != qh)
  {
  for (r2=-1; r2<=1; r2++)
    for (c2=-1; c2<=1; c2++)
      {
          if (r2 == 0  &&  c2 == 0)
              continue;
          newR=queue[qt]/128+r2;
          newC=queue[qt]%128+c2;
          if ((newR) <=0  ||  (newR) >= 128-distance  || (newC) <= 0  ||  (newC) >= 128-distance)
              continue;
          if (seed_matrix[(newR)*128+newC]!=0)
              continue;
		  //test criteria to join region
		  diffAngle=getAngle(averageX,averageY,averageZ,Sx[newR*128+newC],Sy[newR*128+newC],Sz[newR*128+newC]);
          if (abs(diffAngle) > angle_threshold)
              continue;
          //otherwise join pixel to current region...
          count++;
          averageX=(averageX*(count-1)+Sx[newR*128+newC])/count;
          averageY=(averageY*(count-1)+Sy[newR*128+newC])/count;
          averageZ=(averageZ*(count-1)+Sz[newR*128+newC])/count;
          seed_matrix[newR*128+newC]=1;
          CopyImage[newR*128+newC]=nextLevelIntensity;
          queue[qh]=newR*128+newC;
          qh=(qh+1)%MAX_QUEUE;
      }
  qt=(qt+1)%MAX_QUEUE;
  }
  nextLevelIntensity+=50;
}

void convertTo3D()
{
int	r,c,ROWS=128,COLS=128;
double	cp[7];
double	xangle,yangle,dist;
double	SlantCorrection;

cp[0]=1220.7;		/* horizontal mirror angular velocity in rpm */
cp[1]=32.0;		/* scan time per single pixel in microseconds */
cp[2]=(COLS/2)-0.5;		/* middle value of columns */
cp[3]=1220.7/192.0;	/* vertical mirror angular velocity in rpm */
cp[4]=6.14;		/* scan time (with retrace) per line in milliseconds */
cp[5]=(ROWS/2)-0.5;		/* middle value of rows */
cp[6]=10.0;		/* standoff distance in range units (3.66cm per r.u.) */

cp[0]=cp[0]*3.1415927/30.0;	/* convert rpm to rad/sec */
cp[3]=cp[3]*3.1415927/30.0;	/* convert rpm to rad/sec */
cp[0]=2.0*cp[0];		/* beam ang. vel. is twice mirror ang. vel. */
cp[3]=2.0*cp[3];		/* beam ang. vel. is twice mirror ang. vel. */
cp[1]/=1000000.0;		/* units are microseconds : 10^-6 */
cp[4]/=1000.0;			/* units are milliseconds : 10^-3 */

  for (r=0; r<ROWS; r++)
    {
    for (c=0; c<COLS; c++)
      {
      SlantCorrection=cp[3]*cp[1]*((double)c-cp[2]);
      xangle=cp[0]*cp[1]*((double)c-cp[2]);
      yangle=(cp[3]*cp[4]*(cp[5]-(double)r))+SlantCorrection;
      dist=(double)RangeImage[r*COLS+c]+cp[6];
      Pz[r*COLS+c]=sqrt((dist*dist)/(1.0+(tan(xangle)*tan(xangle))+(tan(yangle)*tan(yangle))));
      Px[r*COLS+c]=tan(xangle)*Pz[r*COLS+c];
      Py[r*COLS+c]=tan(yangle)*Pz[r*COLS+c];
      }
    }
}
