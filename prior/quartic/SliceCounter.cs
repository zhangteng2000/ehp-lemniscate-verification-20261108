using System;
using System.IO;
using System.Collections.Generic;
using System.Threading.Tasks;

// Integer proof kernel. All real input enclosures and model comparisons are made by slice_data.py.
public class SliceCounter
{
    const int Mask = (1 << 25) - 1;
    static int nx, directions, depth, nodes;
    static int[] data;
    static ParallelOptions pool = new ParallelOptions();
    class State { public byte[] bound; public bool full; }
    static Dictionary<int,State> cache = new Dictionary<int,State>();
    class Shape
    {
        public bool full;
        public long[][] gen;
        public long[] rad, crude, clip, clipRad;
        public long sourceError, clipError;
    }
    struct Range
    {
        public int lo, hi;
        public Range(int l, int h) { lo=l; hi=h; }
    }
    static void Controls(long[] t, int sample, int node, int level)
    {
        long max=0;
        for(int i=0;i<81;i++) max=Math.Max(max,Math.Abs(t[i]));
        if (max >= (1L<<60)) throw new Exception("source guard");
        int shift=8;
        while((max>>shift)>=Mask) shift++;
        int off=(sample*nodes+node)*81;
        for(int j=0;j<9;j++)
            for(int i=0;i<9;i++) data[off+9*j+i]=(int)(t[9*i+j]>>shift);
        if(level==depth)return;
        long[] left=new long[81], right=new long[81];
        long[] a=new long[9];
        for(int i=0;i<9;i++)
        {
            Array.Copy(t,9*i,a,0,9);
            left[9*i]=a[0];right[9*i+8]=a[8];
            for(int k=1;k<=8;k++)
            {
                for(int j=0;j<=8-k;j++) a[j]=(a[j]+a[j+1])>>1;
                left[9*i+k]=a[0]; right[9*i+8-k]=a[8-k];
            }
        }
        Controls(left,sample,2*node+1,level+1);
        Controls(right,sample,2*node+2,level+1);
    }
    static long Dot(long[] g, int p)
    {
        return g[0]*data[p]+g[1]*data[p+1]+g[2]*data[p+2]+g[3]*data[p+3]+
               g[4]*data[p+4]+g[5]*data[p+5]+g[6]*data[p+6]+g[7]*data[p+7];
    }
    static long AbsDot(long[] g, int p)
    {
        return g[0]*Math.Abs(data[p])+g[1]*Math.Abs(data[p+1])+g[2]*Math.Abs(data[p+2])+g[3]*Math.Abs(data[p+3])+
               g[4]*Math.Abs(data[p+4])+g[5]*Math.Abs(data[p+5])+g[6]*Math.Abs(data[p+6])+g[7]*Math.Abs(data[p+7]);
    }
    static int Sign(Shape s,int p)
    {
        long m=Dot(s.gen[0],p), am=Math.Abs(m);
        if(am<=s.sourceError)return 0;
        if(am>s.sourceError+AbsDot(s.crude,p))return m>0?1:-1;
        long err=s.sourceError+AbsDot(s.rad,p);
        for(int i=1;i<s.gen.Length && err<am;i++)err+=Math.Abs(Dot(s.gen[i],p));
        return am>err?(m>0?1:-1):0;
    }
    // 1 all inside the clipping window, -1 all outside, 0 unresolved.
    static int Clip(Shape s,int p)
    {
        bool inside=true, outside=true;
        for(int j=0;j<9 && (inside || outside);j++,p+=9)
        {
            long val=s.clip[4]*data[p+8],err=s.clipError;
            for(int k=0;k<4;k++){val+=s.clip[k]*data[p+k];err+=s.clipRad[k]*Math.Abs(data[p+k]);}
            err+=s.clipRad[4]*Math.Abs(data[p+8]);
            if(val+err>=0) inside=false;
            if(val-err<=0) outside=false;
        }
        return inside?1:outside?-1:0;
    }
    static Range Count(Shape s,int sample,int node,int level,bool inside)
    {
        int start=(sample*nodes+node)*81;
        int prev=0,loss=0,first=0,last=0;
        for(int j=0;j<9;j++)
        {
            int v=Sign(s,start+9*j);
            if(j==0)first=v;
            if(j==8)last=v;
            if(v!=0)
            {
                if((j&1)!=0) v=-v;
                if(prev!=0 && prev!=v)loss++;
                prev=v;
            }
        }
        int hi=8-loss;
        if(hi==0)return new Range(0,0);
        if(!inside)
        {
            int c=Clip(s,start);
            if(c<0)return new Range(0,0);
            inside=c>0;
        }
        int lo=(inside && first*last<0)?1:0;
        if(level==depth || hi==lo)return new Range(lo,hi);
        Range a=Count(s,sample,2*node+1,level+1,inside);
        Range b=Count(s,sample,2*node+2,level+1,inside);
        return new Range(Math.Max(lo,a.lo+b.lo),Math.Min(hi,a.hi+b.hi));
    }
    static void Run(string line)
    {
        string[] parts=line.Split(' ');
        if(parts[0]=="D"){cache.Remove(int.Parse(parts[1]));Console.WriteLine("OK");return;}
        int at=0;
        int id=int.Parse(parts[at++]),parent=int.Parse(parts[at++]);
        Shape s=new Shape(); s.full=parts[at++]=="1";
        int n=int.Parse(parts[at++]);s.gen=new long[n][];
        for(int l=0;l<n;l++)
        {
            s.gen[l]=new long[8];
            for(int a=0;a<8;a++)s.gen[l][a]=long.Parse(parts[at++]);
        }
        s.rad=new long[8];s.crude=new long[8];
        long norm=0;
        for(int a=0;a<8;a++)
        {
            s.crude[a]=s.rad[a]=long.Parse(parts[at++]);
            if(s.rad[a]<0)throw new Exception("remainder sign");
            for(int l=1;l<n;l++)s.crude[a]+=Math.Abs(s.gen[l][a]);
            norm+=Math.Abs(s.gen[0][a])+s.crude[a];
        }
        if(norm<=0 || norm> (1L<<60)/(Mask+3))throw new Exception("coefficient guard");
        s.sourceError=2*norm;
        s.clip=new long[5];s.clipRad=new long[5];norm=0;
        for(int a=0;a<5;a++)s.clip[a]=long.Parse(parts[at++]);
        for(int a=0;a<5;a++)
        {
            s.clipRad[a]=long.Parse(parts[at++]);
            if(s.clipRad[a]<0)throw new Exception("clip error");
            norm+=Math.Abs(s.clip[a])+s.clipRad[a];
        }
        if(parts.Length!=at || norm> (1L<<60)/(Mask+3))throw new Exception("clip guard");
        s.clipError=2*norm;
        State old=null;
        if(parent>=0 && cache.ContainsKey(parent) && cache[parent].full==s.full)old=cache[parent];
        State state=new State{bound=new byte[nx*directions],full=s.full};
        int[] hiSum=new int[directions],loSum=new int[directions];
        Parallel.For(0,directions*16,pool,(int piece)=>
        {
            int dr=piece/16,start=(piece%16)*nx/16,end=(piece%16+1)*nx/16;
            int hs=0,ls=0;
            for(int x=start;x<end;x++)
            {
                int sample=dr*nx+x;
                int prior=old==null?8:old.bound[sample];
                int lo=prior>>4,hi=prior&15;
                if(lo!=hi)
                {
                    Range r=Count(s,sample,0,0,s.full);
                    lo=Math.Max(lo,r.lo);hi=Math.Min(hi,r.hi);
                    if(s.full){lo+=lo&1; hi-=hi&1;}
                    if(hi<lo)throw new Exception("inconsistent slice count");
                }
                state.bound[sample]=(byte)((lo<<4)+hi);
                hs+=hi;ls+=lo;
            }
            System.Threading.Interlocked.Add(ref hiSum[dr],hs);
            System.Threading.Interlocked.Add(ref loSum[dr],ls);
        });
        cache[id]=state;
        Console.WriteLine(String.Join(",",hiSum)+" "+String.Join(",",loSum));
    }
    public static void Main(string[] args)
    {
        depth=int.Parse(args[1]);nodes=(1<<(depth+1))-1;
        if(depth>10 || depth<0)throw new Exception("depth");
        pool.MaxDegreeOfParallelism=int.Parse(args[2]);
        using(var b=new BinaryReader(File.OpenRead(args[0])))
        {
            nx=b.ReadInt32();directions=b.ReadInt32();
            if(b.ReadInt32()!=36 || b.ReadInt32()!=4)throw new Exception("table format");
            data=new int[checked(nx*directions*nodes*81)];
            for(int sample=0;sample<nx*directions;sample++)
            {
                long[] t=new long[81];
                for(int j=0;j<81;j++) t[j]=b.ReadInt64();
                Controls(t,sample,0,0);
            }
            if(b.BaseStream.Position!=b.BaseStream.Length)throw new Exception("table length");
        }
        Console.WriteLine("READY");
        string line;
        while((line=Console.ReadLine())!=null)Run(line);
    }
}

