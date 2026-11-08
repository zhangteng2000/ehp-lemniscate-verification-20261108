using System;
using System.IO;
using System.Globalization;
using System.Collections.Generic;
using System.Threading.Tasks;

// Every arithmetic operation that could overflow is checked at compile time.
// Bounds are integers in units 2^-24; right shifts on signed long round down.
public class JointCubic {
    const int Bits=24;
    const long S=1L<<Bits;
    struct IV {
        public long L,H;
        public IV(long l,long h) {L=l;H=h;if(l>h)throw new Exception("bad interval");}
        public static IV Int(long v){return new IV(v*S,v*S);}
        public static IV Point(long v){return new IV(v,v);}
        public static IV operator +(IV a, IV b){return new IV(a.L+b.L,a.H+b.H);}
        public static IV operator -(IV a, IV b){return new IV(a.L-b.H,a.H-b.L);}
        public static IV operator *(IV a, IV b){
            long p=a.L*b.L,q=a.L*b.H,r=a.H*b.L,s=a.H*b.H;
            long lo=Math.Min(Math.Min(p,q),Math.Min(r,s)),hi=Math.Max(Math.Max(p,q),Math.Max(r,s));
            return new IV(lo>>Bits,-((-hi)>>Bits));
        }
        public IV Times(long c){return c>=0?new IV(L*c,H*c):new IV(H*c,L*c);}
        public IV Div(long d){
            if(d<=0)throw new Exception();
            return new IV(L/d-(L%d<0?1:0),H/d+(H%d>0?1:0));
        }
        public IV Sq(){
            long lo=Math.Max(0,Math.Max(L,-H)),hi=Math.Max(-L,H);
            return new IV(lo*lo>>Bits,-((-hi*hi)>>Bits));
        }
    }
    class Cell {
        public int Q,A,B,Count;
        public Cell(int q,int a,int b){Q=q;A=a;B=b;}
    }
    static IV[] sin;
    static int N,Qn,An,Bn,K,M,words;
    static bool first;
    static long top,Hmax,hs;
    static IV Sin(int j){return sin[j%(2*N)];}
    static IV Cos(int j){return Sin(j+N/2);}
    static int Pop(ulong a){
        unchecked {
            a-= (a>>1)&0x5555555555555555UL;
            a=(a&0x3333333333333333UL)+((a>>2)&0x3333333333333333UL);
            return (int)((((a+(a>>4))&0x0f0f0f0f0f0f0f0fUL)*0x0101010101010101UL)>>56);
        }
    }
    static long[] Errors(int iq){
        // Use the largest shape parameter in any incident cell.
        IV v=IV.Point(top*Math.Min(Qn,iq+1)/Qn),f=first?IV.Int(1):v,g=first?v:IV.Int(1);
        IV f2=f.Sq(),g2=g.Sq(),w=g2+(f2*f).Times(4).Div(27);
        IV step=IV.Point(top/Qn);
        long[] err=new long[M+1];
        for(int j=0;j<=M;j++){
            long h0=Math.Abs(-Hmax+j*hs);
            IV r=IV.Point(Math.Min(Hmax,h0+hs)),small=IV.Point(Math.Max(0,h0-hs));
            IV r2=r.Sq(),r3=r2*r,r4=r2.Sq();
            IV ang=w.Times(484)+(f*g*r).Times(338)+(f2*r2).Times(96)+(g*r3).Times(484)+(f*r4).Times(16);
            IV alpha=w.Times(16)+(f*g*r).Times(8)+(f2*r2).Times(16).Div(3)+(g*r3).Times(16);
            IV beta=g2.Times(4)+(f*g*r).Times(2)+(g*r3).Times(4);
            IV shape=first?IV.Int(2):(r2.Times(2)+f.Times(8).Div(9));
            IV hess=(f*r2).Times(48)+(g*r).Times(24)+f2.Times(2)-small.Sq().Sq().Times(120);
            hess=IV.Point(Math.Max(0,hess.H));
            //  pi^2 < 10. Coordinates theta,alpha,beta cover pi,pi/2,pi.
            IV total=ang.Times(10).Div(8L*K*K)+alpha.Times(10).Div(32L*An*An)
                +beta.Times(10).Div(8L*Bn*Bn)+(shape*step.Sq()).Div(8)
                +(hess*IV.Point(hs).Sq()).Div(8);
            err[j]=total.H;
        }
        return err;
    }
    static ulong[] Vertex(int iq,int ia,int ib,long[] err){
        long raw=top*iq/Qn;
        IV f=first?IV.Int(1):IV.Point(raw),g=first?IV.Point(raw):IV.Int(1);
        int alpha=ia*(N/(2*An)),beta=ib*(N/Bn);
        ulong[] output=new ulong[K*words],prev=new ulong[words],cur=new ulong[words];
        for(int k=0;k<=K;k++){
            int t=k*(N/K);
            IV X=Cos(7*t+alpha),Y=Sin(7*t+alpha),b=f*Sin(5*t+alpha),c=g*Sin(4*t+alpha+beta);
            IV zi=f*Sin(2*t),yc=Y*c;
            IV h4=(f*Cos(2*t)).Times(4),h3=(X*c).Times(4),h2=zi.Sq()-b.Sq().Times(4).Div(3),
               h1=(yc*zi).Times(2),h0=yc.Sq()+(Y*b*b*b).Times(4).Div(27);
            Array.Clear(cur,0,words);bool last=false;
            for(int j=0;j<=M;j++){
                long zh=-Hmax+j*hs,h=Math.Abs(zh),s3=zh<0?-h3.H:h3.L,s1=zh<0?-h1.H:h1.L;
                long v=((-4*h*h)>>Bits)+h4.L;
                v=((v*h)>>Bits)+s3;
                v=((v*h)>>Bits)+h2.L;
                v=((v*h)>>Bits)+s1;
                v=((v*h)>>Bits)+h0.L;
                bool ok=v>err[j];
                if(ok && last)cur[(j-1)>>6]|=1UL<<((j-1)&63);
                last=ok;
            }
            if(k>0)for(int w=0;w<words;w++)output[(k-1)*words+w]=prev[w]&cur[w];
            ulong[] tmp=prev;prev=cur;cur=tmp;
        }
        return output;
    }
    static ulong[] Fetch(Dictionary<int,ulong[]> memo,int iq,int ia,int ib,int dq,long[] err){
        int key=dq*(Bn+1)+ib;ulong[] r;
        if(!memo.TryGetValue(key,out r)){r=Vertex(iq,ia,ib,err);memo.Add(key,r);}
        return r;
    }
    static void Slab(List<Cell> cells){
        if(cells.Count==0)return;
        int q=cells[0].Q,lastA=-2;
        long[][] errors={Errors(q),Errors(q+1)};
        Dictionary<int,ulong[]> lo=new Dictionary<int,ulong[]>(),hi=new Dictionary<int,ulong[]>();
        ulong[][] v=new ulong[8][];
        foreach(Cell c in cells){
            if(c.A!=lastA){
                if(c.A==lastA+1)lo=hi;else lo=new Dictionary<int,ulong[]>();
                hi=new Dictionary<int,ulong[]>();lastA=c.A;
            }
            for(int dq=0;dq<2;dq++)for(int da=0;da<2;da++)for(int db=0;db<2;db++)
                v[dq*4+da*2+db]=Fetch(da==0?lo:hi,q+dq,c.A+da,c.B+db,dq,errors[dq]);
            int count=0;
            for(int k=0;k<K*words;k++)
                count+=Pop(v[0][k]&v[1][k]&v[2][k]&v[3][k]&v[4][k]&v[5][k]&v[6][k]&v[7][k]);
            c.Count=count;
        }
    }
    public static void Main(string[] args){
        if(args.Length!=10)throw new Exception("kind Q A B K M input output trig workers");
        first=args[0]=="A";if(!first && args[0]!="B")throw new Exception("kind");
        Qn=int.Parse(args[1]);An=int.Parse(args[2]);Bn=int.Parse(args[3]);
        K=int.Parse(args[4]);M=int.Parse(args[5]);
        top=first?5*S/16:9*S/4;Hmax=first?5*S/4:7*S/4;words=(M+63)/64;
        if(top%Qn!=0 || (2*Hmax)%M!=0)throw new Exception("non-dyadic mesh");
        hs=2*Hmax/M;
        string[] lines=File.ReadAllLines(args[8]);string[] header=lines[0].Split(' ');
        N=int.Parse(header[0]);
        if(int.Parse(header[1])!=Bits||N%K!=0||N%(2*An)!=0||N%Bn!=0||lines.Length!=2*N+1)
            throw new Exception("angle grid");
        sin=new IV[2*N];for(int j=0;j<sin.Length;j++){
            string[] r=lines[j+1].Split(' ');sin[j]=new IV(long.Parse(r[0]),long.Parse(r[1]));
        }
        var byQ=new List<Cell>[Qn];for(int q=0;q<Qn;q++)byQ[q]=new List<Cell>();
        if(args[6]=="all"){
            for(int q=0;q<Qn;q++)for(int a=0;a<An;a++)for(int b=0;b<Bn;b++)byQ[q].Add(new Cell(q,a,b));
        }else foreach(string line in File.ReadLines(args[6])){
            string[] r=line.Split(' ');int q=int.Parse(r[0]),a=int.Parse(r[1]),b=int.Parse(r[2]);
            if(q<0||q>=Qn||a<0||a>=An||b<0||b>=Bn)throw new Exception("parameter outside");
            byQ[q].Add(new Cell(q,a,b));
        }
        foreach(var slab in byQ)slab.Sort((a,b)=>a.A!=b.A?a.A.CompareTo(b.A):a.B.CompareTo(b.B));
        Parallel.For(0,Qn,new ParallelOptions{MaxDegreeOfParallelism=int.Parse(args[9])},q=>Slab(byQ[q]));
        int solved=0,used=0,min=int.MaxValue;
        using(var wr=new StreamWriter(args[7])){
            foreach(var slab in byQ)foreach(var c in slab){
                used++;if(first?3330L*c.Count>742L*K*M:1998L*c.Count>530L*K*M)solved++;
                min=Math.Min(min,c.Count);wr.WriteLine(c.Q+" "+c.A+" "+c.B+" "+c.Count);
            }
        }
        Console.WriteLine(args[0]+" "+Qn+"x"+An+"x"+Bn+" K="+K+" M="+M+" proved="+solved+"/"+used+" min_cells="+min);
    }
}
