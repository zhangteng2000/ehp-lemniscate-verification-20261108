using System;
using System.IO;
using System.Collections.Generic;
using System.Threading.Tasks;

// A second evaluation of the five-dimensional covering. A target leaf
// is proved by its own area test, irrespective of the producer's count.
public class SexticCubicAudit {
    const double Pad=1.0/33554432.0;
    static double[] sine,hs,hs2;
    static int Q,A,B,T,H,W;
    static bool chartA;
    static double end,radius;
    class Leaf { public int q,a,b,count; }
    static double Sn(int k){return sine[k&8191];}
    static double Cs(int k){return Sn(k+2048);}
    static double[] Budget(int q) {
        double vmax=(end/Q)*Math.Min(Q,q+1),f=chartA?1:vmax,g=chartA?vmax:1;
        double dt=1.0/T,da=0.5/A,db=1.0/B,ds=end/Q,dh=2*radius/H;
        double[] limit=new double[H+1];
        for(int i=0;i<=H;i++) {
            double hi=Math.Min(radius,Math.Abs(hs[i])+dh),lo=Math.Max(0,Math.Abs(hs[i])-dh);
            double f2=f*f,r2=hi*hi,r3=r2*hi,r4=r2*r2,w=g*g+4.0/27*f2*f;
            double angle=324*w+242*f*g*hi+64*f2*r2+324*g*r3+16*f*r4;
            double alpha=16*w+8*f*g*hi+16.0/3*f2*r2+16*g*r3;
            double beta=4*g*g+2*f*g*hi+4*g*r3;
            double shape=chartA?2:2*r2+8.0/9*f;
            double height=Math.Max(0,48*f*r2+24*g*hi+2*f2-120*(lo*lo)*(lo*lo));
            // 79/8 > (355/113)^2 > pi^2.
            limit[i]=((79.0/8)*(dt*dt*angle+da*da*alpha+db*db*beta)
                     +shape*ds*ds+height*dh*dh)/8+Pad;
        }
        return limit;
    }
    static ulong[] Mask(int q,int a,int b,double[] lim) {
        double f=chartA?1:(end/Q)*q,g=chartA?(end/Q)*q:1;
        int phase=a*(2048/A),offset=b*(4096/B);
        ulong[] result=new ulong[T*W],prev=new ulong[W];
        for(int row=0;row<=T;row++) {
            int v=row*(4096/T);
            double y=Sn(6*v+phase),x=Cs(6*v+phase),bi=f*Sn(4*v+phase);
            double c=g*Sn(3*v+phase+offset),zi=f*Sn(2*v),zr=f*Cs(2*v);
            double yc=y*c,free=(4.0/27)*y*bi*bi*bi,odd=4*x*c,even=(4.0/3)*bi*bi;
            bool last=false;
            for(int m=0;m<W;m++) {
                ulong bits=0;
                // Assign a face to its right height endpoint (0,1] ... (H-1,H].
                if(m==0){
                    double h=hs[0],h2=hs2[0],s=h*zi+yc;
                    last=s*s+free+h2*(h*odd-even+4*h2*(zr-h2))>lim[0];
                }
                for(int n=0;n<64;n++) {
                    int j=1+64*m+n;
                    double h=hs[j],h2=hs2[j],s=h*zi+yc;
                    double det=s*s+free+h2*(h*odd-even+4*h2*(zr-h2));
                    bool good=det>lim[j];
                    if(last&&good)bits|=1UL<<n;
                    last=good;
                }
                if(row!=0)result[(row-1)*W+m]=bits&prev[m];
                prev[m]=bits;
            }
        }
        return result;
    }
    static void Slab(List<Leaf> leaves) {
        if(leaves.Count==0)return;
        int q=leaves[0].q;double[][] lim={Budget(q),Budget(q+1)};
        var cache=new Dictionary<long,ulong[]>();int old=-3;
        foreach(Leaf l in leaves) {
            if(l.a!=old){
                var keep=new Dictionary<long,ulong[]>();
                foreach(var item in cache)if(item.Key/(2*(B+1))==l.a)keep.Add(item.Key,item.Value);
                cache=keep;old=l.a;
            }
            var sets=new ulong[8][];
            for(int a=0;a<2;a++)for(int b=0;b<2;b++)for(int j=0;j<2;j++){
                long key=(2*(l.a+a)+j)*(B+1)+l.b+b;
                ulong[] v;
                if(!cache.TryGetValue(key,out v)){
                    v=Mask(q+j,l.a+a,l.b+b,lim[j]);cache.Add(key,v);
                }
                sets[4*a+2*b+j]=v;
            }
            int total=0;
            for(int i=0;i<sets[0].Length;i++){
                ulong mask=sets[0][i]&sets[1][i]&sets[2][i]&sets[3][i]&sets[4][i]&sets[5][i]&sets[6][i]&sets[7][i];
                for(int k=0;k<8;k++){total+=pop[(int)(mask&255)];mask>>=8;}
            }
            l.count=total;
        }
    }
    static byte[] pop=new byte[256];
    public static void Main(string[] args) {
        if(args.Length!=10)throw new Exception("chart Q A B theta height passed-leaves audit-output sine-table workers");
        chartA=args[0]=="A";if(!chartA&&args[0]!="B")throw new Exception();
        Q=int.Parse(args[1]);A=int.Parse(args[2]);B=int.Parse(args[3]);T=int.Parse(args[4]);H=int.Parse(args[5]);W=H/64;
        if(H%64!=0||4096%T!=0||2048%A!=0||4096%B!=0)throw new Exception();
        for(int i=1;i<256;i++)pop[i]=(byte)(pop[i/2]+(i&1));
        end=chartA?5.0/16:2.25;radius=chartA?1.25:1.75;
        var tr=File.ReadAllLines(args[8]);if(tr.Length!=8193||tr[0]!="4096 46")throw new Exception();
        sine=new double[8192];for(int i=0;i<8192;i++){
            var z=tr[i+1].Split(' ');long x=long.Parse(z[0]),y=long.Parse(z[1]);
            if(y<x||y-x>1)throw new Exception();
            sine[i]=x/70368744177664.0;
        }
        hs=new double[H+1];hs2=new double[H+1];
        for(int j=0;j<=H;j++){hs[j]=(2*j-H)*(radius/H);hs2[j]=hs[j]*hs[j];}
        var byQ=new List<Leaf>[Q];for(int q=0;q<Q;q++)byQ[q]=new List<Leaf>();
        foreach(string line in File.ReadLines(args[6])){
            var z=line.Split(' ');var l=new Leaf{q=int.Parse(z[0]),a=int.Parse(z[1]),b=int.Parse(z[2])};
            if(l.q<0||l.q>=Q||l.a<0||l.a>=A||l.b<0||l.b>=B)throw new Exception();
            byQ[l.q].Add(l);
        }
        foreach(var ls in byQ)ls.Sort((k,l)=>k.a!=l.a?k.a.CompareTo(l.a):k.b.CompareTo(l.b));
        Parallel.For(0,Q,new ParallelOptions{MaxDegreeOfParallelism=int.Parse(args[9])},q=>Slab(byQ[q]));
        using(var f=new StreamWriter(args[7]))foreach(var ls in byQ)foreach(var l in ls)
            f.WriteLine(l.q+" "+l.a+" "+l.b+" "+l.count);
    }
}
