using System;
using System.IO;
using System.Collections.Generic;
using System.Threading.Tasks;

// Independent verifier of integer slice-count certificates.
public sealed class CheckSlices
{
    static int N,K,D,B;
    static int[] net;
    static ParallelOptions cfg = new ParallelOptions();
    sealed class Frame { public byte[] rows; public bool full; }
    static Dictionary<int,Frame> frames = new Dictionary<int,Frame>();
    sealed class Bounds
    {
        public long[][] g;
        public long[] rad,cm,cr;
        public long eps,ceps;
        public bool full;
    }
    static void Sub(long[,] v,int row,int index,int level)
    {
        long mx=0;
        foreach(long x in v) mx=Math.Max(mx,Math.Abs(x));
        if(mx>=(1L<<60))throw new InvalidDataException("source size");
        int power=8;
        while((mx>>power)>=(1<<25)-1)++power;
        int start=(row*B+index)*81;
        for(int k=0;k<9;k++)for(int j=0;j<9;j++)
            net[start+k*9+j]=checked((int)(v[k,j]>>power)+1);
        if(level==D)return;
        long[,] left=new long[9,9],right=new long[9,9];
        for(int k=0;k<9;k++)
        {
            long[] values=new long[9];
            for(int j=0;j<9;j++)values[j]=v[k,j];
            int width=8;
            for(int step=0;step<=8;step++)
            {
                left[k,step]=values[0];right[k,8-step]=values[width];
                if(width==0)break;
                for(int j=0;j<width;j++)values[j]=(values[j]+values[j+1])>>1;
                width--;
            }
        }
        Sub(left,row,2*index+1,level+1);Sub(right,row,2*index+2,level+1);
    }
    static int Allowed(Bounds b,int start,int j)
    {
        long center=0,delta=b.eps;
        for(int i=0;i<8;i++)
        {
            long p=net[start+9*i+j];
            center+=b.g[0][i]*p;delta+=b.rad[i]*Math.Abs(p);
        }
        for(int t=1;t<b.g.Length;t++)
        {
            long value=0;
            for(int i=0;i<8;i++)value+=b.g[t][i]*net[start+9*i+j];
            delta+=Math.Abs(value);
        }
        return center>delta?1:center < -delta?-1:0;
    }
    static int Clip(Bounds b,int start)
    {
        bool positive=true,negative=true;
        for(int j=0;j<9 && (positive||negative);j++)
        {
            long c=0,err=b.ceps;
            for(int h=0;h<5;h++)
            {
                long p=net[start+(h==4?72:9*h)+j];c+=b.cm[h]*p;err+=b.cr[h]*Math.Abs(p);
            }
            positive &= c>err;negative &= c < -err;
        }
        return positive?1:negative?-1:0;
    }
    // Low count in high nibble, high count in low nibble. The upper comes
    // from the two-state sign recurrence, not the forced-alternation formula.
    static byte Box(Bounds b,int row,int node,int level,bool inside)
    {
        int start=(row*B+node)*81,first=Allowed(b,start,0),last=first;
        int pos=first<0?-32:0,neg=first>0?-32:0;
        for(int j=1;j<=8;j++)
        {
            last=Allowed(b,start,j);
            int a=last<0?-32:Math.Max(pos,neg+1),c=last>0?-32:Math.Max(neg,pos+1);
            pos=a;neg=c;
        }
        int upper=Math.Max(pos,neg);
        if(upper==0)return 0;
        if(!inside)
        {
            int c=Clip(b,start);
            if(c>0)return 0;
            inside=c<0;
        }
        int lower=inside && first*last<0?1:0;
        if(level<D && upper>lower)
        {
            int r=Box(b,row,2*node+2,level+1,inside),l=Box(b,row,2*node+1,level+1,inside);
            upper=Math.Min(upper,(r&15)+(l&15));lower=Math.Max(lower,(r>>4)+(l>>4));
        }
        if(upper<lower)throw new InvalidDataException("root interval");
        return (byte)(16*lower+upper);
    }
    static void Request(string str)
    {
        string[] t=str.Split(' ');
        if(t[0]=="D"){frames.Remove(int.Parse(t[1]));Console.WriteLine("OK");return;}
        int pos=0,id=int.Parse(t[pos++]),par=int.Parse(t[pos++]);
        Bounds b=new Bounds();b.full=int.Parse(t[pos++])==1;
        int ng=int.Parse(t[pos++]);b.g=new long[ng][];
        for(int i=0;i<ng;i++){
            b.g[i]=new long[8];for(int j=0;j<8;j++){long v=long.Parse(t[pos++]);b.g[i][j]=v;b.eps+=Math.Abs(v);}
        }
        b.rad=new long[8];
        for(int j=0;j<8;j++){long v=long.Parse(t[pos++]);if(v<0)throw new InvalidDataException();b.rad[j]=v;b.eps+=v;}
        b.cm=new long[5];b.cr=new long[5];
        for(int j=0;j<5;j++){b.cm[j]=long.Parse(t[pos++]);b.ceps+=Math.Abs(b.cm[j]);}
        for(int j=0;j<5;j++){b.cr[j]=long.Parse(t[pos++]);if(b.cr[j]<0)throw new InvalidDataException();b.ceps+=b.cr[j];}
        if(pos!=t.Length||Math.Max(b.eps,b.ceps)>((1L<<60)/((1<<25)+3)))throw new InvalidDataException("integer range");
        Frame old;
        if(!frames.TryGetValue(par,out old)||old.full!=b.full)old=null;
        Frame current=new Frame{ rows=new byte[N*K],full=b.full };
        int[] hi=new int[K],lo=new int[K];
        Parallel.For(0,K*32,cfg,(int task)=>{
            int k=task/32,s=task%32,a=0,h=0;
            for(int j=N*s/32;j<N*(s+1)/32;j++)
            {
                int row=k*N+j,prev=old==null?8:old.rows[row],l=prev>>4,u=prev&15;
                if(u!=l){
                    int r=Box(b,row,0,0,b.full);
                    u=Math.Min(u,r&15);l=Math.Max(l,r>>4);
                    if(b.full){u=2*(u/2);l=2*((l+1)/2);}
                    if(u<l)throw new InvalidDataException("inherited count");
                }
                current.rows[row]=(byte)(16*l+u);a+=l;h+=u;
            }
            System.Threading.Interlocked.Add(ref hi[k],h);
            System.Threading.Interlocked.Add(ref lo[k],a);
        });
        frames[id]=current;
        Console.WriteLine(String.Join(",",hi)+" "+String.Join(",",lo));
    }
    public static void Main(string[] args)
    {
        D=int.Parse(args[1]);cfg.MaxDegreeOfParallelism=int.Parse(args[2]);
        if(D<0||D>10)throw new InvalidDataException();
        B=(1<<(D+1))-1;
        using(BinaryReader r=new BinaryReader(File.OpenRead(args[0])))
        {
            N=r.ReadInt32();K=r.ReadInt32();if(r.ReadInt32()!=36||r.ReadInt32()!=4)throw new InvalidDataException();
            net=new int[checked(N*K*B*81)];
            for(int row=0;row<N*K;row++)
            {
                long[,] a=new long[9,9];
                for(int k=0;k<9;k++)for(int j=0;j<9;j++)a[k,j]=r.ReadInt64();
                Sub(a,row,0,0);
            }
            if(r.BaseStream.Position!=r.BaseStream.Length)throw new InvalidDataException("trailing source data");
        }
        Console.WriteLine("READY");
        string line;while((line=Console.ReadLine())!=null)Request(line);
    }
}

