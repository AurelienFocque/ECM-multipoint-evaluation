#include <pari/pari.h>
#include <pari/paripriv.h>

// Returns product tree of x-v[i] mod n
// n doesn't need to be prime
GEN
FpX_product_tree(GEN v,GEN n)
{
	int nbetages=0;
	int t=1;
	int d=lg(v)-2;
	int j;
	int prop1;
	while(d>0)
	{
		d=d/2;
		nbetages++;
		t=2*t;
	}
	GEN tmp=cgetg(nbetages+1,t_VEC);
	GEN g0=cgetg(t+1,t_VEC);

	prop1=(t/(t-lg(v)+1));
	if(prop1*(t-lg(v)+1)!=t){prop1++;}
	int k=0;
	int l=0;
	while(l<=t-1)
	{   
		if(k<=lg(v)-2&&l%prop1!=0)
		{    
			gel(g0,1+l)=mkpoln(2,gen_1,gneg(gel(v,1+k)));
			k++;
			l++;
		}
		else
		{
			gel(g0,1+l)=mkpoln(2,gen_0,mkintn(1,1));
			l++;
		}


	}

	gel(tmp,1)=g0;
	j=1;
	GEN g1;
	while(j<nbetages)
	{
		t=t/2;
		g1=cgetg(t+1,t_VEC);
		k=0;
		while(k<=t-1)
		{
			gel(g1,1+k)=FpX_mul(gel(gel(tmp,j),1+2*k),gel(gel(tmp,j),2+2*k),n);
			k++;
		}
		gel(tmp,1+j)=g1;
		j++;
	}
	return tmp;

}
//return the multipoint evaluation of P modulo n.
//NB : impossible divisions with FpX_rem will never occur because all the members of the product tree are monic polynomials,
//so n doesn't need to be prime
GEN
FpX_multipoint_eval(GEN P,GEN listepoints,GEN n)
{   
	GEN res;
	GEN arbre=FpX_product_tree(listepoints,n);
	pari_printf("Arbre des sous-produits : %Ps\n",arbre);
	int k=lg(arbre)-1;
	GEN tmp1,tmp2;

	tmp1=cgetg(3,t_VEC);
	gel(tmp1,1)=FpX_rem(P,gel(gel(arbre,k),1),n);
	gel(tmp1,2)=FpX_rem(P,gel(gel(arbre,k),2),n);
	k--;
	int j =4;
	int i;
	while(k>0)
	{
		tmp2=tmp1;
		tmp1=cgetg(j+1,t_VEC);
		i=1;
		while(i<=j)
		{
			gel(tmp1,i)=FpX_rem(gel(tmp2,(i-1)/2+1),gel(gel(arbre,k),i),n);
			i++;
		}
	k--;
	j=2*j;
	}

	//Now we clear the artificial zeros
	int size=lg(listepoints)-1;
	res=cgetg(size+1,t_VEC);
	int l=1;
	int m =1;
	while(l<=size)
	{
		if(RgX_equal(gel(gel(arbre,1),m),mkpoln(1,gen_1))==0)
		{
		//cast them into scalars
			if(lg(gel(tmp1,m))<3)//beware of zero polynomial
				gel(res,l)=gen_0;
			else
				gel(res,l)=gel(gel(tmp1,m),2);
			l++;

		}
		m++;
	}
	gerepileupto(avma,res);
	return res;
}
GEN
ECM_Fast_product(GEN a,GEN b,GEN n)
// ECM stage 2 requires de compute prod(1<=i<j<=B)(a_i-b_j) mod N and then take its gcd with n
//a and b stand for the list of coeffients a_i and b_i
{   
	GEN res=gen_1;
	GEN list=FpX_product_tree(a,n);
	int k = lg(list);
	GEN F=FpX_mul(gel(gel(list,k-1),1),gel(gel(list,k-1),2),n);
	GEN tmp=FpX_multipoint_eval(F,b,n);
	k=lg(tmp)-1;
	int i =1;
	while(i<=k)
	{
		res=Fp_mul(res,gel(tmp,i),n);
		i++;
	}
	gerepileupto(avma,res);
	return res;
}

int
main()
{

  pari_init(10000000000,2);
  GEN P=RgX_shift(mkpoln(1,gen_1),50);
  GEN n=mkintn(1,101);
  GEN v=cgetg(102,t_VEC);
  int k=1;
  while(k<=101)
  {
    gel(v,k)=mkintn(1,k);
    k++;
  }
  pari_printf("Exemple (symbole de Legendre multipoint):\nÉvaluation multipoint de %Ps en les points %Ps modulo %Ps\n",P,v,n);
  pari_printf("Évaluations : %Ps\n",FpX_multipoint_eval(P,v,n));
  gerepileupto(avma,0);
  pari_close();
  return 0;
}
