event plot (i >= 0)
{
    if (dimension == 1 && matlab_out)
    {        
        char name[11];

        /* Write header file */
        sprintf(name,"fort.t%04d",Frame);
        FILE *fp = fopen(name,"w");
        fprintf(fp,"%20.16f %20s\n",t,"time");
        fprintf(fp,"%10d %20s\n",meqn,"meqn");
        fprintf(fp,"%10d %20s\n",1,"ngrids");
        fprintf(fp,"%10d %20s\n",maux,"maux");
        fprintf(fp,"%10d %20s\n",dimension,"dim");
        fclose(fp);

        /* Write data file */
        sprintf(name,"fort.q%04d",Frame);
        fp = fopen(name,"w");
        fprintf(fp,"%10d %20s\n",1,"grid_number");
        fprintf(fp,"%10d %20s\n",1,"AMR_level");
        fprintf(fp,"%10d %20s\n",N,"mx");
        fprintf(fp,"%24.16f %20s\n",X0,"xlow");  /* xlow */

        /* Assume uniform Cartesian grid */
        double dx = L0/N;
        fprintf(fp,"%24.16f %20s\n",dx,"dx");
        fprintf(fp,"\n");
        foreach()
        {
            for(scalar s in statevars)
                fprintf(fp,"%24.16e",s[]);
            fprintf(fp,"\n");
        }
        fclose(fp);        
    }
    if (dimension == 2 && matlab_out)
    {        
        char name[11];

        /* Write header file */
        sprintf(name,"fort.t%04d",Frame);
        FILE *fp = fopen(name,"w");
        fprintf(fp,"%20.16f %20s\n",t,"time");
        fprintf(fp,"%10d %20s\n",meqn,"meqn");
        fprintf(fp,"%10d %20s\n",1,"ngrids");
        fprintf(fp,"%10d %20s\n",maux,"maux");
        fprintf(fp,"%10d %20s\n",dimension,"dim");
        fclose(fp);

        /* Write out face and vertex data for Matlab 'patch' */
        FILE *fpv, *fpf, *fpc;

        sprintf(name,"fort.v%04d",Frame);
        fpv = fopen(name,"wb");

        sprintf(name,"fort.f%04d",Frame);
        fpf = fopen(name,"wb");

        sprintf(name,"fort.c%04d",Frame);
        fpc = fopen(name,"wb");

        double xpadd[5] = {0,1,1,0,0};
        double ypadd[5] = {0,0,1,1,0};
        int n = 0;
        foreach()
        {                
            double h = Delta;
            double xlow = x-h/2.0;
            double ylow = y-h/2.0;
            double vp[2*5];
            int f[5];
            for(int i = 0; i < 5; i++)
            {                
                f[i]        = i + 5*n + 1;
                vp[2*i]     = xlow + h*xpadd[i];
                vp[2*i + 1] = ylow + h*ypadd[i];
            }
            fwrite(vp,sizeof(double),10,fpv);
            fwrite(f,sizeof(int),5,fpf);

            int m = 0;
            double cdata[meqn];
            for(scalar s in statevars)
            {
                cdata[m] = s[];
                m++;
            }
            fwrite(cdata,sizeof(double),meqn,fpc);
            n++;
        }
        fclose(fpc);
        fclose(fpv);
        fclose(fpf);
    }    
    printf("Output Frame %d at time t = %12.4e\n",Frame, t);
    printf("\n");
    Frame++;
}
