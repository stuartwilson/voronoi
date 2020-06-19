#include "voronoi.h"
using morph::Tools;

int main(int argc, char **argv){

    if(argc<3){ cout<<"Supply number of points and number of temperatures. Exiting..."<<endl; return 0; }

    // Construct Voronoi diagram
    Voronoi V(stoi(argv[1]));

    // Get vector of Voronoi cells
    vector<cell> Q;
    for(int i=0;i<V.N;i++){ Q.push_back(cell(V,i)); }

    double MinSegLen = 1e9;
    double MaxSegLen = -1;
    for(int i=0;i<Q.size();i++){
        for(int j=0;j<Q[i].segL.size();j++){
            if(MaxSegLen<Q[i].segL[j]){
                MaxSegLen = Q[i].segL[j];
            }
            if(MinSegLen>Q[i].segL[j]){
                MinSegLen = Q[i].segL[j];
            }
        }
    }

    cout<<"Max: "<<MaxSegLen<<endl;
    cout<<"Min: "<<MinSegLen<<endl;

    // Gibb's sampling (following MacKay, 2003, p401)


    int I = 100000;     // iterations
    int N = V.N;
    vector<int> S;      // state of each cell

    double kB = 1.0;    // Bolzmann constant

    double J = 1.0; // coupling
    double H = 0.0; // applied field

    int nTemp = stoi(argv[2]);

    vector<double> Tout(nTemp);
    vector<double> Mout(nTemp);
    vector<double> Eout(nTemp);
    vector<double> Vout(nTemp);
    vector<double> Sout(nTemp);

    double maxTemp = 50.0;
    double minTemp = 0.1;
    double tempIncrement = (maxTemp-minTemp)/((double)nTemp-1);
    vector<double> Temps(nTemp);
    for(int i=0;i<nTemp;i++){
        Temps[i] = maxTemp - tempIncrement*(double)i;
    }

    // reset (no reason not to do this on every test?? But following MacKay for now)
    S.resize(N,-1);
    for(int i=0;i<N;i++){
        if(Tools::randDouble()<0.5){
            S[i] = 1;
        }
    }


    for(int t=0;t<nTemp;t++){
        double T = Temps[t];
        double beta = 2./(kB * T); // note 2 as $S\in\pm 1$ (1 for $S\in 1,0$)
        // Ising dynamics
        double M = 0.;
        int count = 0;
        vector<double> E;
        for(int i=0;i<I;i++){
            int k = floor(Tools::randDouble()*N); // choose random spin
            //double F = 0.; for(int j=0;j<Q[k].N;j++){ F += J*S[Q[k].G[j]] + H; } // field
            double F = 0.;
            for(int j=0;j<Q[k].N;j++){
                double J = Q[k].segLscaled[j];
                F += J*S[Q[k].G[j]] + H;
            } // field
            double P = 1./(1+exp(-beta*F)); // prob +1
            if(Tools::randDouble()<P){ S[k] = 1; } else { S[k] = -1; } // update state

            if(i>I/3){

                // Magnetization
                double m = 0.;
                for(int l=0;l<N;l++){ m+=S[l]; }
                m /= (double)N;
                M += m;

                // Energy
                double f = 0.;
                for(int l=0;l<N;l++){
                    //for(int j=0;j<Q[l].N;j++){ f += 0.5*J*S[l]*S[Q[l].G[j]] + H; } // field
                    for(int j=0;j<Q[l].N;j++){
                       double J = Q[k].segLscaled[j];
                       //double J = 1.0;
                        f += J*S[l]*S[Q[l].G[j]] + H;
                    } // field
                }

                f = -f*0.5/double(N);
                E.push_back(f);

                count++;
            }
        }

        // Standard deviation of energy
        double meanE = 0.;
        for(int q=0;q<E.size();q++){
            meanE += E[q];
        }
        meanE /= (double)E.size();

        double stdE = 0.;
        for(int q=0;q<E.size();q++){
            stdE += (E[q]-meanE)*(E[q]-meanE);
        }
        stdE = pow(stdE/((double)E.size()-1.0),0.5);

        // Store values
        Tout[t] = T;
        Mout[t] = M/(double)count;
        Eout[t] = meanE;
        Vout[t] = stdE;
        cout<<"T: "<<Tout[t]<<", Energy: "<<Eout[t]<<endl;
    }

    HdfData outdata("ising.h5");
    outdata.add_contained_vals ("T", Tout);
    outdata.add_contained_vals ("M", Mout);
    outdata.add_contained_vals ("E", Eout);
    outdata.add_contained_vals ("V", Vout);




    // save output to recreate plot
    if(true){
        V.save_points();
        Q[5].save();
    }

}

// VARYING TESSELATION SIZE SEEMS TO WORK FOR LOWER N, BUT IT COULD BE JUST THAT AT HIGHER N THERE IS AN INCREASING CHANCE OF HAVING SOME PATHALOGICAL FEATURE IN THE TESSELATION THAT OTHERWISE COULD BE ACCOUNTED FOR TO ENABLE HIGHER N SIMS??
