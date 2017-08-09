#include "adjust.h"

void ADJUST(int Attemp, int Nacc, double &Dmax, int Attv, int Accv, double &Vmax, double Succ, double HBOX[2]){

        static int attempp = 0;
        static int naccp   = 0;
        static int attvp   = 0;
        static int accvp   = 0;
        double dro, vmaxo, frac;
//	printf("%s , HBOX: %6.3f Dmax: %6.3f, Attemp: %6d, Nacc: %6d\n", __func__, HBOX, Dmax, Attemp, Nacc);
//---displacement:
        if ( Attemp ==0 || attempp >= Attemp ){
                naccp   = Nacc;
                attempp = Attemp;
        }else{
                frac = double(Nacc-naccp)/double(Attemp-attempp);
                dro  = Dmax;
                Dmax = Dmax*fabs( frac/(Succ) );
//		printf("%s , dro: %6.3f Dmax: %6.3f, frac: %6.3f\n", __func__, dro, Dmax, frac);
//              ---limit the change:
                if( Dmax/dro > 1.5 ) Dmax = dro*1.5 ;
                if( Dmax/dro < 0.5 ) Dmax = dro*0.5 ;
                if( Dmax > HBOX[0]/2.0) Dmax = HBOX[0]/2.0;
//                printf("Max. displ. set to: %6.3f, ( old : %6.3f )\n", Dmax, dro);
//                printf("Frac. acc. : %6.3f, attempts : %6d, succes:  %6d \n", frac, (Attemp-attempp), (Nacc-naccp));
//              ---store nacc and attemp for next use           
                naccp   = Nacc;
                attempp = Attemp;
        }

//---volume:
        if ( Attv ==0 || attvp >= Attv ){
                accvp = Accv;
                attvp = Attv;
        }else{
                frac = double(Accv-accvp)/double(Attv-attvp);
                vmaxo = Vmax;
                Vmax = Vmax*fabs( frac/(Succ) );
//              ---limit the change:
                if( Vmax/vmaxo > 1.5 ) Vmax = vmaxo*1.5 ;               
                if( Vmax/vmaxo < 1.5 ) Vmax = vmaxo*0.5 ;
//                printf("Max. vol. chan. set to: %6.3f, ( old : %6.3f )\n", Vmax, vmaxo); 
//                printf("Frac. acc. : %6.3f, attempts : %6d, succes:  %6d \n", frac, (Attv-attvp), (Accv-accvp) );
                accvp = Accv;
                attvp = Attv;
        } 

}

