import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.io.*;
import java.util.*;
import java.security.*;
import javax.swing.*;
import javax.imageio.ImageIO;

// ------------- class Point ------------------------------
class Pnt {
    public int x,y;
    public Pnt() {};
    public Pnt(int x1, int y1) {
        x = x1;
        y = y1;
    }
    public boolean equals(Pnt other) {
        return (x == other.x && y == other.y);
    }
}

// ------------- class CutTheRoots itself --------------
public class CutTheRootsVis {
    final int SZ = 1024;            // field size
    int NP, NR, AVGR;               // number of plants, number of roots, average roots per plant
    ArrayList<Pnt> p = new ArrayList<Pnt>(); // coordinates of points (fixed)
    ArrayList<Integer> rlen = new ArrayList<Integer>(); // Length of parent root
    ArrayList<Double> angle = new ArrayList<Double>(); // Angle of root
    ArrayList<Integer> plant = new ArrayList<Integer>(); // Root belongs to this plant
    int[] dead;                     // List of dead points
    int[] roots;                    // List of edges for the roots
    int[] pointsPar;                // Point parameter
    int[] cX1, cY1, cX2, cY2;       // Cuts
    double rootsLength;             // Length of roots for all plants
    double[] rootDead;
    // ---------------------------------------------------
    void generate(String seed) {
      try {
        SecureRandom rnd = SecureRandom.getInstance("SHA1PRNG");
        rnd.setSeed(Long.parseLong(seed));
        // generate points by sampling each coordinate uniformly, without duplicates
        int i, j, k;
        // number of points
        if (seed.equals("1"))
        {
            NP = 5;
            AVGR = 100;
            NR = NP*AVGR;
        } else {
            NP = rnd.nextInt(100) + 5;
            AVGR = 10 + rnd.nextInt(1000);
            NR = NP*AVGR;
        }
        System.out.println("NP = " + NP + " NR = " + NR + " AVGR = " + AVGR);
        roots = new int[NR*2];
        rootsLength = 0;
        // generate the plants
        boolean ok;
        for (i = 0; i < NP; ++i)
        {
            do {
                Pnt apnt = new Pnt(rnd.nextInt(SZ), rnd.nextInt(SZ));
                ok = true;
                if ( (apnt.x-SZ/2)*(apnt.x-SZ/2) + (apnt.y-SZ/2)*(apnt.y-SZ/2) > SZ*SZ/4 ) ok = false;
                for (j = 0; j < i && ok; ++j)
                    if (apnt.equals(p.get(j)))
                        ok = false; 
                if (ok)
                {
                    plant.add(i);
                    p.add(apnt);
                    rlen.add(rnd.nextInt(100) + 10);
                    angle.add(-100.0);
                }
            }
            while (!ok);
        }
        // generate the roots
        for (i = 0; i < NR; ++i)
        {
            int r = rnd.nextInt(p.size());
            int l = 1 + (rlen.get(r) * (50+rnd.nextInt(40)) / 100);
            double scale = 1.0;
            do {
                int dx,dy;
                double ang = angle.get(r);
                if (ang<-50)
                {
                    ang = rnd.nextDouble()*3.1415926*2.0;
                } else
                {
                    ang += scale*(-1.0 + rnd.nextDouble()*2.0);
                }
                dx = p.get(r).x + (int)(Math.cos(ang)*l);
                dy = p.get(r).y + (int)(Math.sin(ang)*l);
                ok = false;
                if ( (dx-SZ/2)*(dx-SZ/2) + (dy-SZ/2)*(dy-SZ/2) > SZ*SZ/4 )
                {
                    l /= 2;
                    scale *= 2.0;
                    continue;
                }
                roots[i*2] = r;
                roots[i*2+1] = p.size();
                rootsLength += Math.sqrt( (p.get(r).x - dx)*(p.get(r).x - dx) +
                                          (p.get(r).y - dy)*(p.get(r).y - dy) );
                p.add(new Pnt(dx, dy));
                rlen.add(l);
                angle.add(ang);
                plant.add(plant.get(r));
                ok = true;
            }
            while (!ok);
        }

        System.out.println("Length of all roots = " + rootsLength);

        // convert points to parameter array
        pointsPar = new int[2 * p.size()];
        for (i = 0; i < p.size(); ++i) {
            pointsPar[2 * i] = p.get(i).x;
            pointsPar[2 * i + 1] = p.get(i).y;
        }

      }
      catch (Exception e) { 
        addFatalError("An exception occurred while generating test case.");
        e.printStackTrace(); 
      }
    }
    // ---------------------------------------------------
    int lineIntersection(double Ax, double Ay, double Bx, double By,
                         double Cx, double Cy, double Dx, double Dy,
                         double[] XY)
    {
        double  distAB, theCos, theSin, tmp, ABpos;
        //  Fail if either line is undefined.
        if ((Math.abs(Ax-Bx)<1e-10 && Math.abs(Ay-By)<1e-10) ||
            (Math.abs(Cx-Dx)<1e-10 && Math.abs(Cy-Dy)<1e-10)) return 0;

        //  (1) Translate the system so that point A is on the origin.
        Bx-=Ax; By-=Ay;
        Cx-=Ax; Cy-=Ay;
        Dx-=Ax; Dy-=Ay;

        //  Discover the length of segment A-B.
        distAB = Math.sqrt(Bx*Bx+By*By);

        //  (2) Rotate the system so that point B is on the positive X axis.
        theCos = Bx/distAB;
        theSin = By/distAB;
        tmp  = Cx*theCos+Cy*theSin;
        Cy   = Cy*theCos-Cx*theSin;
        Cx   = tmp;
        tmp  = Dx*theCos+Dy*theSin;
        Dy   = Dy*theCos-Dx*theSin;
        Dx   = tmp;

        //  (3) Discover the position of the intersection point along line A-B.
        ABpos = Dx+(Cx-Dx)*Dy/(Dy-Cy);

        //  Fail if segment C-D doesn't cross line A-B.
        if ((Cy<-1e-10 && Dy<-1e-10) || (Cy>1e-10 && Dy>1e-10)) return 0;

        //  (3) Discover the position of the intersection point along line A-B.
        ABpos = Dx+(Cx-Dx)*Dy/(Dy-Cy);

        //  Fail if segment C-D crosses line A-B outside of segment A-B.
        if (ABpos<-1e-10 || ABpos>distAB+1e-10) return 0;

        //  (4) Apply the discovered position to line A-B in the original coordinate system.
        XY[0]=Ax+ABpos*theCos;
        XY[1]=Ay+ABpos*theSin;
        XY[2]=ABpos/distAB;
        return 1;
    }
    // ---------------------------------------------------
    double calcScore() {

        double rootLen = 0;
        double[] XY = new double[3];

        rootDead = new double[NR];
        for (int i=0;i<NR;i++) rootDead[i] =0;
        dead = new int[p.size()];
        for (int i=0;i<p.size();i++) dead[i] =0;
        for (int i=0;i<NR;i++)
        {
            // two points of this root
            int p1 = roots[i*2];
            int p2 = roots[i*2+1];
            if (dead[p1]==1)
            {
                rootDead[i] = 1.0;
                // this root is dead, don't add it to the length
                dead[p1] = dead[p2] = 1; // propagate the effect
                continue;
            }
            // check if any cut cuts this root
            double[] cutXY = new double[3];
            cutXY[2] = 1e20;
            for (int j=0;j<cX1.length;j++)
            {
                int res = lineIntersection(p.get(p1).x, p.get(p1).y, p.get(p2).x, p.get(p2).y,
                                 cX1[j], cY1[j], cX2[j], cY2[j], XY);
                if (res==1)
                {
                    if (XY[2]<cutXY[2])
                    {
                        cutXY[2] = XY[2];
                    }
                }
            }
            double totalLen = Math.sqrt( (p.get(p1).x - p.get(p2).x)*(p.get(p1).x - p.get(p2).x) +
                                         (p.get(p1).y - p.get(p2).y)*(p.get(p1).y - p.get(p2).y) );
            if (cutXY[2]<1e10)
            {
                dead[p2] = 1;
                rootLen += cutXY[2] * totalLen;
                rootDead[i] = cutXY[2];
            } else
            {
                // not cut, add the length
                rootLen += totalLen;
            }
        }

        // check that all plants have been separated
        // yes, this is slow, find a faster way yourself, it is part of the challenge ;)
        for (int i1=0;i1<NP;i1++)
            for (int i2=i1+1;i2<NP;i2++)
            {
                boolean ok = false;
                for (int c=0;c<cX1.length;c++)
                {
                    int res = lineIntersection(p.get(i1).x, p.get(i1).y, p.get(i2).x, p.get(i2).y,
                                                 cX1[c], cY1[c], cX2[c], cY2[c], XY);
                    if (res==1)
                    {
                        ok = true;
                        break;
                    }
                }
                if (!ok)
                {
                    addFatalError("Plants "+i1+" and "+i2+" not separated with a cut.");
                    return 0;
                }
            }

        // now, if all are valid, score is always non-0
        System.out.println("Length of alive roots = " + rootLen);
        double score = 1000000.0 * rootLen / rootsLength;
        return score;
    }
    // ---------------------------------------------------
    public double runTest(String seed) {
      try {
        int i, j;
        generate(seed);
        double score = 0;
        if (proc != null) {
            // get the return and parse it
            int[] ret;
            try { ret = makeCuts(NP, pointsPar, roots); }
            catch (Exception e) {
                addFatalError("Failed to get result from makeCuts.");
                return 0;
            }

            int ncuts = ret.length/4;
            if (ret.length%4!=0)
            {
                addFatalError("Return should contain a multiple of four elements.");
                return 0;
            }
            if (ncuts>NP*4)
            {
                addFatalError("Maximum number of cuts allowed is "+(NP*4)+", you returned "+ncuts+" cuts.");
                return 0;
            }
            for (i=0;i<ret.length;i++)
            {
                if (ret[i]<0 || ret[i]>SZ)
                {
                    addFatalError("Return values should in the range of [0,"+SZ+"]. Your value was "+ret[i]+".");
                    return 0;
                }
            }
            cX1 = new int[ncuts];
            cY1 = new int[ncuts];
            cX2 = new int[ncuts];
            cY2 = new int[ncuts];
            int idx = 0;
            for (i=0;i<ncuts;i++)
            {
                cX1[i] = ret[idx++];
                cY1[i] = ret[idx++];
                cX2[i] = ret[idx++];
                cY2[i] = ret[idx++];
                if (cX1[i]==cX2[i] && cY1[i]==cY2[i])
                {
                    addFatalError("Cut points be at least 1 unit apart.");
                    return 0;
                }
                // extend lines to start and end outside region
                int dx = (cX2[i]-cX1[i]);
                int dy = (cY2[i]-cY1[i]);
                if (dx==0)
                {
                    cY1[i] = 0;
                    cY2[i] = SZ;
                } else if (dy==0)
                {
                    cX1[i] = 0;
                    cX2[i] = SZ;
                } else
                {
                    if (cX1[i]>cX2[i])
                    {
                        int tmp = cX1[i];
                        cX1[i] = cX2[i];
                        cX2[i] = tmp;
                        tmp = cY1[i];
                        cY1[i] = cY2[i];
                        cY2[i] = tmp;
                        dx = -dx;
                        dy = -dy;
                    }
                    int ox = cX1[i];
                    int oy = cY1[i];
                    int extendL = 1+Math.max(ox/dx, oy/dy);
                    int extendR = 1+Math.max((SZ-ox)/dx, (SZ-oy)/dy);
                    cX1[i] = ox - dx*extendL;
                    cY1[i] = oy - dy*extendL;
                    cX2[i] = ox + dx*extendR;
                    cY2[i] = oy + dy*extendR;
                }
            }
            score = calcScore();
        }

        if (vis) {
            // draw the image
            jf.setSize(SZX,SZY);
            jf.setVisible(true);
            draw();
        }
        if (saveFile!=null)
            saveCase(saveFile);

        return score;
      }
      catch (Exception e) { 
        addFatalError("An exception occurred while trying to process your program's results.");
        e.printStackTrace(); 
        return 0.0;
      }
    }
// ------------- visualization part ----------------------
    static String exec;
    static boolean vis;
    static Process proc;
    static String saveFile = null;
    JFrame jf;
    Vis v;
    InputStream is;
    OutputStream os;
    BufferedReader br;
    // problem-specific drawing params
    final int SZX = SZ+2,SZY=SZ+2;
    volatile boolean ready;
    volatile int Ncur;
    volatile int[] Pcur;
    int[][] coordToPoint;
    // ---------------------------------------------------
    int[] makeCuts(int NP_, int[] points, int[] roots) throws IOException
    {   // pass the params to the solution and get the return
        int i;
        StringBuffer sb = new StringBuffer();
        sb.append(NP_).append('\n');
        sb.append(points.length).append('\n');
        for (i = 0; i < points.length; ++i)
            sb.append(points[i]).append('\n');
        sb.append(roots.length).append('\n');
        for (i = 0; i < roots.length; ++i)
            sb.append(roots[i]).append('\n');
        os.write(sb.toString().getBytes());
        os.flush();
        // get the return - an array of strings
        int nret = Integer.parseInt(br.readLine());
        //System.out.println(nret);
        int[] ret = new int[nret];
        for (i = 0; i < nret; ++i)
            ret[i] = Integer.parseInt(br.readLine());
        return ret;
    }
    // ---------------------------------------------------
    void draw() {
        if (!vis) return;
        v.repaint();
    }

    BufferedImage drawCase()
    {
        int i;
        BufferedImage bi = new BufferedImage(SZX+10,SZY+10,BufferedImage.TYPE_INT_RGB);
        Graphics2D g2 = (Graphics2D)bi.getGraphics();
        //background
        g2.setColor(new Color(0xD3D3D3));
        g2.fillRect(0,0,SZX+10,SZY+10);
        g2.setColor(Color.WHITE);
        g2.fillRect(0,0,SZ+1,SZ+1);
        //frame
        g2.setColor(Color.BLACK);
        g2.drawRect(0,0,SZ+1,SZ+1);

        g2.setColor(Color.GRAY);
        g2.fillOval(0,0,SZ,SZ);
        if (p.size()>NP)
        {
            for (i=0;i<NR;i++)
            {
                float hue = (float)(plant.get(roots[i*2])) / NP;
                int p1 = roots[i*2];
                int p2 = roots[i*2+1];
                if (rootDead==null || dead[p2]==0)
                {
                    g2.setColor(Color.getHSBColor(hue, 0.9f, 1.0f));
                    g2.drawLine(p.get(roots[i*2]).x, p.get(roots[i*2]).y, p.get(roots[i*2+1]).x, p.get(roots[i*2+1]).y);
                } else
                {
                    if (dead[p1]==0)
                    {
                        g2.setColor(Color.getHSBColor(hue, 0.9f, 1.0f));
                        int nx = p.get(roots[i*2]).x+(int)((p.get(roots[i*2+1]).x-p.get(roots[i*2]).x)*rootDead[i]);
                        int ny = p.get(roots[i*2]).y+(int)((p.get(roots[i*2+1]).y-p.get(roots[i*2]).y)*rootDead[i]);
                        g2.drawLine(p.get(roots[i*2]).x, p.get(roots[i*2]).y, nx, ny);
                        g2.setColor(Color.getHSBColor(hue, 0.9f, 0.2f));
                        g2.drawLine(nx, ny, p.get(roots[i*2+1]).x, p.get(roots[i*2+1]).y);
                    } else
                    {
                        g2.setColor(Color.getHSBColor(hue, 0.9f, 0.2f));
                        g2.drawLine(p.get(roots[i*2]).x, p.get(roots[i*2]).y, p.get(roots[i*2+1]).x, p.get(roots[i*2+1]).y);
                    }
                }
            }
            for (i=0;i<NP;i++)
            {
                float hue = (float)(i) / NP;
                g2.setColor(Color.getHSBColor(hue, 0.9f, 1.0f));
                g2.fillOval(p.get(i).x-5, p.get(i).y-5, 10, 10);
            }
            g2.setColor(Color.WHITE);
            for (i=0;i<NP;i++)
            {
                g2.fillOval(p.get(i).x-2, p.get(i).y-2, 4, 4);
            }
        }

        // draw cuts
        if (cX1!=null)
        for (i=0;i<cX1.length;i++)
        {
            g2.setColor(Color.WHITE);
            g2.drawLine((int)cX1[i], (int)cY1[i], (int)cX2[i], (int)cY2[i]);
        }
        return bi;
    }

    void saveCase(String filename)
    {
        try {
            BufferedImage bi = drawCase();
            File outputfile = new File(filename);
            ImageIO.write(bi, "png", outputfile);
        } catch (IOException e) {
            addFatalError("An exception occurred while trying to save the image.");
            e.printStackTrace();
        }
    }
    // ---------------------------------------------------
    public class Vis extends JPanel implements WindowListener {
        public void paint(Graphics g) {
          try {
            BufferedImage bi = drawCase();
            g.drawImage(bi,0,0,SZX+10,SZY+10,null);
      }
      catch (Exception e) { e.printStackTrace(); }
        }
        public Vis() {
            jf.addWindowListener(this);
        }
    // ---------------------------------------------------
        //WindowListener
        public void windowClosing(WindowEvent e) { 
            if(proc != null)
                try { proc.destroy(); } 
                catch (Exception ex) { ex.printStackTrace(); }
            System.exit(0); 
        }
        public void windowActivated(WindowEvent e) { }
        public void windowDeactivated(WindowEvent e) { }
        public void windowOpened(WindowEvent e) { }
        public void windowClosed(WindowEvent e) { }
        public void windowIconified(WindowEvent e) { }
        public void windowDeiconified(WindowEvent e) { }
    }
    // ---------------------------------------------------
    public CutTheRootsVis(String seed) {
        //interface for runTest
        if (vis)
        {   jf = new JFrame();
            v = new Vis();
            jf.getContentPane().add(v);
        }
        if (exec != null) {
            try {
                Runtime rt = Runtime.getRuntime();
                proc = rt.exec(exec);
                os = proc.getOutputStream();
                is = proc.getInputStream();
                br = new BufferedReader(new InputStreamReader(is));
                new ErrorReader(proc.getErrorStream()).start();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        System.out.println("Score = "+runTest(seed));
        if (proc != null)
            try { proc.destroy(); } 
            catch (Exception e) { e.printStackTrace(); }
    }
    // ---------------------------------------------------
    public static void main(String[] args) {
        String seed = "1";
        vis = false;
        for (int i = 0; i<args.length; i++)
        {   if (args[i].equals("-seed"))
                seed = args[++i];
            if (args[i].equals("-exec"))
                exec = args[++i];
            if (args[i].equals("-save"))
                saveFile = args[++i];
            if (args[i].equals("-vis"))
                vis = true;
        }
        CutTheRootsVis f = new CutTheRootsVis(seed);
    }
    // ---------------------------------------------------
    void addFatalError(String message) {
        System.out.println(message);
    }
}

class ErrorReader extends Thread{
    InputStream error;
    public ErrorReader(InputStream is) {
        error = is;
    }
    public void run() {
        try {
            byte[] ch = new byte[50000];
            int read;
            while ((read = error.read(ch)) > 0)
            {   String s = new String(ch,0,read);
                System.out.print(s);
                System.out.flush();
            }
        } catch(Exception e) { }
    }
}
