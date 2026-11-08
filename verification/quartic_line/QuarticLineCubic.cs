using System;
using System.IO;
using System.Threading.Tasks;
using System.Globalization;

// Outward fixed-point arithmetic; compile with /checked+ or -checked+.
// The only unchecked arithmetic is the standard bit population count.
public class QuarticLineCubic {
    const int Bits=24;
    const long Scale=1L<<Bits;
    struct I {
        public long L,H;
        public I(long l,long h){if(l>h)throw new Exception("empty interval");L=l;H=h;}
        public static I Point(long x){return new I(x,x);}
        public static I Int(long x){return Point(x*Scale);}
        public static I operator +(I a,I b){return new I(a.L+b.L,a.H+b.H);}
        public static I operator -(I a,I b){return new I(a.L-b.H,a.H-b.L);}
        public static I operator *(I a,I b){
            long p=a.L*b.L,q=a.L*b.H,r=a.H*b.L,s=a.H*b.H;
            long l=Math.Min(Math.Min(p,q),Math.Min(r,s));
            long h=Math.Max(Math.Max(p,q),Math.Max(r,s));
            return new I(l>>Bits,-((-h)>>Bits));
        }
        public I Times(long x){return x<0?new I(H*x,L*x):new I(L*x,H*x);}
        public I Div(long x){
            if(x<=0)throw new Exception("nonpositive denominator");
            return new I(L/x-(L%x<0?1:0),H/x+(H%x>0?1:0));
        }
        public I Sq(){
            long l=Math.Max(0,Math.Max(L,-H)),h=Math.Max(-L,H);
            return new I((l*l)>>Bits,-(((-h)*h)>>Bits));
        }
    }
    struct C {
        public I Re,Im;
        public C(I a,I b){Re=a;Im=b;}
        public C Rotate(I c,I s){return new C(Re*c-Im*s,Re*s+Im*c);}
    }
    class Input {public int Id,Count;public C E,J,B,Z;public I F,G;}
    static I[] sine;
    static int N,T,H,Words;
    static long Step,Radius=3*Scale;
    static bool factored;
    static I Sin(int j){return sine[j%(2*N)];}
    static I Cos(int j){return Sin(j+N/2);}
    static int Pop(ulong x){unchecked{
        x-=(x>>1)&0x5555555555555555UL;
        x=(x&0x3333333333333333UL)+((x>>2)&0x3333333333333333UL);
        return (int)((((x+(x>>4))&0x0f0f0f0f0f0f0f0fUL)*0x0101010101010101UL)>>56);
    }}
    static long[] Errors(Input a){
        I f=a.F,g=a.G,f2=f.Sq();I W=g.Sq()+(f2*f).Times(4).Div(27);
        long[] err=new long[H+1];
        for(int j=0;j<=H;j++){
            long h=Math.Abs(-Radius+j*Step);
            I r=I.Point(Math.Min(Radius,h+Step)),v=I.Point(Math.Max(0,h-Step));
            I r2=r.Sq(),r3=r2*r,r4=r2.Sq();
            I theta=W.Times(100)+(f*g*r).Times(98)+(f2*r2).Times(56).Div(3)
                +(g*r3).Times(100)+(f*r4).Times(16);
            I hh=(f*r2).Times(48)+(g*r).Times(24)+f2.Times(2)-v.Sq().Sq().Times(120);
            hh=I.Point(Math.Max(0,hh.H));
            // pi^2<10; only theta and height are interpolated.
            I angular=factored?theta.Times(79).Div(64L*T*T):theta.Times(10).Div(8L*T*T);
            err[j]=(angular+(hh*I.Point(Step).Sq()).Div(8)).H;
        }
        return err;
    }
    static void Count(Input a){
        long[] error=Errors(a);ulong[] old=new ulong[Words],cur=new ulong[Words];int count=0;
        bool[] oldDirect=new bool[H+1],curDirect=new bool[H+1];
        for(int k=0;k<=T;k++){
            int t=k*(N/T);
            C e=a.E.Rotate(Cos(4*t),Sin(4*t)),jj=a.J.Rotate(Cos(2*t),Sin(2*t));
            C zz=a.Z.Rotate(Cos(2*t),Sin(2*t));
            I cc=a.B.Re*Sin(t)+a.B.Im*Cos(t),yc=e.Im*cc;
            I h4=zz.Re.Times(4),h3=(e.Re*cc).Times(4),
              h2=zz.Im.Sq()-jj.Im.Sq().Times(4).Div(3),
              h1=(yc*zz.Im).Times(2),h0=yc.Sq()+(e.Im*jj.Im*jj.Im*jj.Im).Times(4).Div(27);
            I free=(e.Im*jj.Im*jj.Im*jj.Im).Times(4).Div(27),even=jj.Im.Sq().Times(4).Div(3);
            Array.Clear(cur,0,Words);bool last=false;
            for(int j=0;j<=H;j++){
                long raw=-Radius+j*Step;
                long lo;
                if(!factored){
                    long h=Math.Abs(raw),c3=raw<0?-h3.H:h3.L,c1=raw<0?-h1.H:h1.L;
                    lo=((-4*h*h)>>Bits)+h4.L;
                    lo=((lo*h)>>Bits)+c3;lo=((lo*h)>>Bits)+h2.L;
                    lo=((lo*h)>>Bits)+c1;lo=((lo*h)>>Bits)+h0.L;
                }else{
                    I h=I.Point(raw),hSquared=h.Sq(),s=h*zz.Im+yc;
                    I d=s.Sq()+free+hSquared*(h*h3-even+(hSquared*(zz.Re-hSquared)).Times(4));
                    lo=d.L;
                }
                bool good=lo>error[j];
                if(factored){
                    curDirect[j]=good;
                    if(j>0&&k>0&&good&&last&&oldDirect[j]&&oldDirect[j-1])count++;
                }else if(good&&last)cur[(j-1)>>6]|=1UL<<((j-1)&63);
                last=good;
            }
            if(!factored&&k>0)for(int j=0;j<Words;j++)count+=Pop(old[j]&cur[j]);
            ulong[] sw=old;old=cur;cur=sw;
            bool[] bs=oldDirect;oldDirect=curDirect;curDirect=bs;
        }
        a.Count=count;
    }
    static C ReadC(string[] a,int i){return new C(new I(long.Parse(a[i]),long.Parse(a[i+1])),new I(long.Parse(a[i+2]),long.Parse(a[i+3])));}
    public static void Main(string[] args){
        if(args.Length!=7)throw new Exception("expanded|factored coefficients output theta height sine workers");
        factored=args[0]=="factored";if(!factored&&args[0]!="expanded")throw new Exception("scheme");
        T=int.Parse(args[3]);H=int.Parse(args[4]);Words=(H+63)/64;
        if(H%64!=0||2*Radius%H!=0)throw new Exception("height mesh");Step=2*Radius/H;
        string[] tr=File.ReadAllLines(args[5]),header=tr[0].Split(' ');N=int.Parse(header[0]);
        int tb=int.Parse(header[1]);
        if(N%T!=0||N%2!=0||tb<Bits||tb>46||tr.Length!=2*N+1)throw new Exception("trig grid");
        sine=new I[2*N];for(int j=0;j<2*N;j++){
            string[] v=tr[j+1].Split(' ');long lo=long.Parse(v[0]),hi=long.Parse(v[1]);
            if(tb>Bits){lo>>=(tb-Bits);hi=-((-hi)>>(tb-Bits));}
            sine[j]=new I(lo,hi);
        }
        string[] source=File.ReadAllLines(args[1]);Input[] data=new Input[source.Length];
        for(int j=0;j<data.Length;j++){
            string[] a=source[j].Split(' ');if(a.Length!=19)throw new Exception("input length");
            data[j]=new Input{Id=int.Parse(a[0]),E=ReadC(a,1),J=ReadC(a,5),B=ReadC(a,9),Z=ReadC(a,13),F=I.Point(long.Parse(a[17])),G=I.Point(long.Parse(a[18]))};
        }
        Parallel.For(0,data.Length,new ParallelOptions{MaxDegreeOfParallelism=int.Parse(args[6])},j=>Count(data[j]));
        using(var w=new StreamWriter(args[2]))foreach(Input a in data)w.WriteLine(a.Id+" "+a.Count);
        Console.WriteLine("Recomputed "+data.Length+" closed quartic line intervals, "+args[0]+" circuit.");
    }
}
