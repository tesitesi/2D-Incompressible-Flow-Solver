/* 
空気力学第二D レポート課題
非圧縮性流体のシミュレーション
03-180327 豊島 拓
*/

import java.util.*;

import javax.lang.model.util.ElementScanner6;

import java.lang.*;
import java.io.*;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.io.IOException;


class Main
{
    private static final int _1 = 1;
	// 初期値設定
	// setflw
	// レイノルズ数
	static final double re = 70.0;
	// クーラン数
	static final double cfl = 0.2;
	// SOR Parameter for p（圧力計算に必要）
	static final double omegap = 1.00;
	static final double maxitp = 100;
	static final double errorp = 1.0e-4;
	// イタレーションの設定
	static final int nlast = 5000;
	static final int nlp = 10;

	// setgrd
	// x方向の計算格子
	static final int mx = 401;
	static final int i1 = 96;
	static final int i2 = 106;
	static final int i3 = 96;
	static final int i4 = 106;
	static final double dx = 1.0 / (double) (i2 - i1);
	// y方向の計算格子
	static final int my = 201;
	static final int j1 = 86;
	static final int j2 = 96;
	static final int j3 = 106;
	static final int j4 = 116;
	static final double dy = 1.0 / (double) (j2 - j1);
	// xy方向の格子を構成
	static final int icent = (i1 + i2) / 2;
	static final int jcent = (j1 + j2) / 2;
	static double[][] x = new double[mx + 2][my + 2];
	static double[][] y = new double[mx + 2][my + 2];

	// slvflw
	static final double dt = cfl * Math.min(dx, dy);
	static double[][] u = new double[mx + 2][my + 2];
	static double[][] v = new double[mx + 2][my + 2];
	static double[][] p = new double[mx + 2][my + 2];
	static double[] cp = new double[my + 2];
	static double time = 0;
	static double cpfore;
	static double cpback;
	static double cpbtm;
	static double cptop;
	static double cd;
	static double cl;
	static double cp1, cp2;
	static int nstep = 0;
	static int nbegin = 0;

	// bcforp
	static int i;
	static int j;

	// poiseq
	static double ux, uy, vx, vy;
	static double[][] rhs = new double[mx + 2][my + 2];
	static double[][] omega = new double[mx + 2][my + 2];
	static double res = 0;
	static double dp = 0;

	// veloeq
	static double[][] urhs = new double[mx + 2][my + 2];
	static double[][] vrhs = new double[mx + 2][my + 2];

	public static void setflw() {
		// 宣言の時点で指定しているためここでは処理しない
	}

	public static void setgrd() {
		for (int i = 1; i <= mx; i++) {
			for (int j = _1; j <= my; j++) {
				x[i][j] = dx * (double) (i - icent);
				y[i][j] = dy * (double)( j - jcent );
	    }
	}
    }
   

    // i+1での u を計算 流れを解く
    public static void slvflw() {
		System.out.print("***    Comp. Conditions:\n");
		System.out.print("CFL     =" + cfl +"\n");
		System.out.print("dt      =" + dt + "\n");
		System.out.print(nlast + " Times Steps to go ...\n");
		System.out.print(">> 2D Incompressible Flow Solver\n");
		System.out.print("Re      =" + re + "\n");
		System.out.print("No. of Grid Points = " + mx + " / " + my + "\n");
		System.out.print("CFL  /   dt    /   Steps \n");
		System.out.print(cfl + " / " + dt + " / " + nlast);
		System.out.println();
		// (1) 初期条件を設定
		intcnd(nbegin,time);
		bcforp();
		bcforv();
		System.out.println("Step	/Res(p) at itr. 	/CD	/CL	/Cp1	/Cp2");
		try {

		FileWriter fws = new FileWriter("./KarmanStat.csv",false);
		PrintWriter pws = new PrintWriter(new BufferedWriter(fws));
		// (2) Time Marching
	for (int n = 1; n <= nlast ; n++) {
	    nstep = n + nbegin;
	    time = time + dt;
	    // <1> 圧力 p に関してポアソン方程式を解く。
	    poiseq();
	    bcforp();
	    // <2> x方向速度 u ,y方向速度 v を更新する。
	    veloeq();
	    bcforv();
	    // <3> 箱周りのCLとCDを計算
	    cd = 0.0;
	    for (int j = j1; j < j2 ; j++) {
		cpfore = ( 2.0 * p[i1][j] + 2.0 * p[i1][j+1] ) / 2.0;
		cpback = ( 2.0 * p[i2][j] + 2.0 * p[i2][j+1] ) / 2.0;
		cd     = cd + ( cpfore - cpback ) * dy;
	    }
	    cl = 0.0;
	    for (int i = i1; i < i2 ; i++) {
		cpbtm = ( 2.0 * p[i][j1] + 2.0 * p[i+1][j1] ) / 2.0;
		cptop = ( 2.0 * p[i][j2] + 2.0 * p[i+1][j2] ) / 2.0;
		cl    = cl + ( cpbtm - cptop ) * dx;
	    }	    // <4> 結果をnlpステップごとにモニター / 出力
	    if (n % nlp == 0) {
		cp1 = 2.0 * p[i2+i2-i1][j1];
		cp2 = 2.0 * p[i2+i2-i1][j2];
		System.out.print(nstep + "	") ;
		System.out.print(res + "	");
		System.out.print(cd + "	");
		System.out.print(cl + "	");
		System.out.print(cp1 + "	");
		System.out.print(cp2);
		System.out.println();
		pws.print(n);
		pws.print(",");
		pws.print(res);
		pws.print(",");
		pws.print(cd);
		pws.print(",");
		pws.print(cl);
		pws.print(",");
		pws.print(cp1);
		pws.print(",");
		pws.print(cp2);
		pws.println();
		}
	}
	pws.close(); 
}
catch (IOException exs) {
		exs.printStackTrace();
	}
// (3) 最終結果の書き出し
	try {
	FileWriter fwp = new FileWriter("./KarmanP.csv", false);
	PrintWriter pwp = new PrintWriter(new BufferedWriter(fwp));
	for (int j = 0; j < my+2; j++) {
	    for (int i = 0; i < mx+2; i++) {
			p[i][j] = 2.0 * p[i][j];
			pwp.print(p[i][j]);
			if (!(i == mx && j == my)) {
				pwp.print(",");
			}
		}
	    }
	System.out.println("Computation Completed.");
	pwp.close();
	} catch (IOException exp) {
	exp.printStackTrace();
	}
}

    // 初期条件の設定
    public static void intcnd(int nbegin, double time) {
    for (int i = 1; i <= mx ; i++) {
	for (int j = 1; j <= my ; j++) {
	    u[i][j] = 1.0;
	    v[i][j] = 0.0;
	    p[i][j] = 0.0;
		//一様流条件を至るところ（壁面上は除く）で与える。
	}
    }
    }
    

    
    // 圧力の境界条件の設定
    public static void bcforp() {
	// (1) 一様流条件（ i = 1, 左端）
	i = 1;
	for (int j = 1; j <= my ; j++) {
	    p[i][j] = 0.0;
	}
	// (2) 下流条件（ i = mx, 右端）
	i = mx;
	for (int j = 1; j <= my ; j++) {
	    p[i][j] = 0.0;
	}
	// (3) 底条件（ j = 1, 下端）
	j = 1;
	for (int i = 1; i <= mx ; i++) {
	    p[i][j] = 0.0;
	}
	// (4) 底条件（ j = my, 上端）
	j = my;
	for (int i = 1; i<= mx; i++) {
	    p[i][j] = 0.0;
	}

	// (5) 壁条件
	// まず四隅を決める -隣の値を流用
	p[i1][j1] = p[i1-1][j1-1];
	p[i1][j2] = p[i1-1][j2+1];
	p[i2][j1] = p[i2+1][j1-1];
	p[i2][j2] = p[i2+1][j2+1];
	i = i1;
	for (int j = j1+1; j < j2 ; j++) {
	    p[i][j] = p[i-1][j];
	}
	i = i2;
	for (int j = j1+1; j < j2 ; j++) {
	    p[i][j] = p[i+1][j];
	}
	j = j1;
	for (int i = i1+1; i < i2 ; i++) {
	    p[i][j] = p[i][j-1];
	}
	j = j2;
	for (int i = i1+1; i < i2 ; i++) {
	    p[i][j] = p[i][j+1];
	}

	p[i3][j3] = p[i3-1][j3-1];
	p[i3][j4] = p[i3-1][j4+1];
	p[i4][j3] = p[i4+1][j3-1];
	p[i4][j4] = p[i4+1][j4+1];
	i = i3;
	for (int j = j3+1; j < j4 ; j++) {
	    p[i][j] = p[i-1][j];
	}
	i = i4;
	for (int j = j3+1; j < j4 ; j++) {
	    p[i][j] = p[i+1][j];
	}
	j = j3;
	for (int i = i3+1; i < i4 ; i++) {
	    p[i][j] = p[i][j-1];
	}
	j = j4;
	for (int i = i3+1; i < i4 ; i++) {
	    p[i][j] = p[i][j+1];
	}
    }
	

    public static void bcforv() {
	// (1) 一様流条件 ( i = 1, 左端 )
	i = 1;
	for (int j = 1; j <= my ; j++) {
	    u[i][j] = 1.0;
	    v[i][j] = 0.0;
	    u[i-1][j] = 1.0;
	    v[i-1][j] = 0.0;
	}
	// (2) 下流条件 ( i = mx, 右端 )
	i = mx;
	for (int j = 1; j <= my ; j++) {
	    u[i][j] = 2.0 * u[i-1][j] - u[i-2][j];
	    v[i][j] = 2.0 * v[i-1][j] - v[i-2][j];
	    u[i+1][j] = 2.0 * u[i][j] - u[i-1][j];
	    v[i+1][j] = 2.0 * v[i][j] - v[i-1][j];
	}
	// 底条件 (j = 1,下端 )
	j = 1;
	for (int i = 1; i <= mx ; i++) {
	    u[i][j] = 2.0 * u[i][j+1] - u[i][j+2];
	    v[i][j] = 2.0 * v[i][j+1] - v[i][j+2];
	    u[i][j-1] = 2.0 * u[i][j] - u[i][j+1];
	    v[i][j-1] = 2.0 * v[i][j] - v[i][j+1];
	}
	// 底条件 (j = my, 上端 )
	j = my;
	for ( int i = 1; i <= mx ; i++) {
	    u[i][j] = 2.0 * u[i][j-1] - u[i][j-2];
	    v[i][j] = 2.0 * v[i][j-1] - v[i][j-2];
	    u[i][j+1] = 2.0 * u[i][j] - u[i][j-1];
	    v[i][j+1] = 2.0 * v[i][j] - v[i][j-1];
	}
	// 壁面条件
	for ( int i = i1; i<= i2 ; i++ ) {
	    for ( int j = j1 ; j<= j2 ; j++ ) {
	    u[i][j] = 0.0;
	    v[i][j] = 0.0;
	    }
	}
	for ( int i = i3; i<= i4 ; i++ ) {
	    for ( int j = j3 ; j<= j4 ; j++ ) {
	    u[i][j] = 0.0;
	    v[i][j] = 0.0;
	    }
	}
	
    }

    // 圧力場を解く
    public static void poiseq() {
	// (1) 右辺の計算
	for (int i = 2; i <= mx-1; i++) {
	    for (int j = 2; j <= my-1; j++) {
		if (outbox(i,j)) { //物体の内部は計算しない
		    ux = ( u[i+1][j] - u[i-1][j] ) / ( 2.0 * dx );
		    uy = ( u[i][j+1] - u[i][j-1] ) / ( 2.0 * dy );
		    vx = ( v[i+1][j] - v[i-1][j] ) / ( 2.0 * dx );
		    vy = ( v[i][j+1] - v[i][j-1] ) / ( 2.0 * dy );
			rhs[i][j] = ( ux + vy ) / dt - ( ux * ux + 2.0 * uy * vx + vy * vy );
			// 渦度もここで計算しておく
			omega[i][j] = vx - uy;
		}
	    }
	}
	// (2) イタレーション
	for (int itr = 1; itr <= maxitp; itr ++) {
	    res = 0.0;
	    for (int i = 2; i <= mx-1; i++) {
		for (int j = 2; j <= my-1; j++) {
		    	if (outbox(i,j)) { //物体の内部は計算しない
			    dp = ( p[i+1][j] + p[i-1][j] ) / ( dx * dx )
				+ ( p[i][j+1] + p[i][j-1] ) / ( dy * dy )
				- rhs[i][j];
			    dp = dp / ( 2.0 / (dx * dx) + 2.0 / (dy * dy) ) - p[i][j];
			    res = res + dp * dp; //残渣
			    p[i][j] = p[i][j] + omegap * dp; // 修正すべき量に緩和係数をかけて修正
			}
		}
	    }
	    bcforp();
	    res = Math.sqrt( res / (double)(mx * my) );
	    if ( res < errorp ) {
		break;
	    }
	}
    }



    // 速度場を解く
    public static void veloeq() {
	// (1) 圧力勾配
        for (int i = 1; i < mx+1; i++) {
	    for (int j = 1; j < my+1; j++) {
	        if (outbox(i,j)) { //物体の内部は計算しない
		    urhs[i][j] = -( p[i+1][j] - p[i-1][j] ) / ( 2.0 * dx );
		    vrhs[i][j] = -( p[i][j+1] - p[i][j-1] ) / ( 2.0 * dy );
		}
	    }
	}

	// (2) 粘性項
        for (int i = 1; i < mx+1; i++) {
	    for (int j = 1; j < my+1; j++) {
	        if (outbox(i,j)) { //物体の内部は計算しない
		    urhs[i][j] = urhs[i][j]
			+ ( u[i+1][j] - 2.0 * u[i][j] + u[i-1][j] ) / ( re * dx * dx )
			+ ( u[i][j+1] - 2.0 * u[i][j] + u[i][j-1] ) / ( re * dy * dy );
		    vrhs[i][j] = vrhs[i][j]
			+ ( v[i+1][j] - 2.0 * v[i][j] + v[i-1][j] ) / ( re * dx * dx )
			+ ( v[i+1][j] - 2.0 * v[i][j] + v[i][j-1] ) / ( re * dy * dy );
		}
	    }
	}

	// (3) x方向のadvenction項
	for (int j = j1+1; j < j2; j++) {
	    u[i1+1][j] = 2.0 * u[i1][j] - u[i1-1][j];
	    u[i2-1][j] = 2.0 * u[i2][j] - u[i2+1][j];
	    v[i1+1][j] = 2.0 * v[i1][j] - v[i1-1][j];
	    v[i2-1][j] = 2.0 * v[i2][j] - v[i2+1][j];
	}
	for (int j = j3+1; j < j4; j++) {
	    u[i3+1][j] = 2.0 * u[i3][j] - u[i3-1][j];
	    u[i4-1][j] = 2.0 * u[i4][j] - u[i4+1][j];
	    v[i3+1][j] = 2.0 * v[i3][j] - v[i3-1][j];
	    v[i4-1][j] = 2.0 * v[i4][j] - v[i4+1][j];
	}
	for (int i = 2; i < mx; i++) {
	    for (int j = 2; j < my; j++) {
	        if (outbox(i,j)) { //物体の内部は計算しない
		    urhs[i][j] = urhs[i][j]
			-u[i][j] * (-u[i+2][j] + 8.0 * ( u[i+1][j] - u[i-1][j] ) + u[i-2][j]) / (12.0 * dx)
			-Math.abs(u[i][j]) * ( u[i+2][j] - 4.0 * u[i+1][j] + 6.0 * u[i][j] - 4.0 * u[i-1][j] + u[i-2][j] ) / ( 4.0 * dx);
		    vrhs[i][j] = vrhs[i][j]
			-u[i][j] * (-v[i+2][j] + 8.0 * ( v[i+1][j] - v[i-1][j] ) + v[i-2][j]) / (12.0 * dx)
			-Math.abs(u[i][j]) * ( v[i+2][j] - 4.0 * v[i+1][j] + 6.0 * v[i][j] - 4.0 * v[i-1][j] + v[i-2][j] ) / ( 4.0 * dx);
		}
	    }
	}

	// (4) y方向のadvection項
	for (int i = i1+1; i < i2; i++) {
	    u[i][j1+1] = 2.0 * u[i][j1] - u[i][j1-1];
	    u[i][j2-1] = 2.0 * u[i][j2] - u[i][j2+1];
	    v[i][j1+1] = 2.0 * v[i][j1] - v[i][j1-1];
	    v[i][j2-1] = 2.0 * v[i][j2] - v[i][j2+1];
	}
	for (int i = i3+1; i < i4; i++) {
	    u[i][j3+1] = 2.0 * u[i][j3] - u[i][j3-1];
	    u[i][j4-1] = 2.0 * u[i][j4] - u[i][j4+1];
	    v[i][j3+1] = 2.0 * v[i][j3] - v[i][j3-1];
	    v[i][j4-1] = 2.0 * v[i][j4] - v[i][j4+1];
	}

	for (int i = 2; i < mx; i++) {
	    for (int j = 2; j < my; j++) {
	        if (outbox(i,j)) { //物体の内部は計算しない
		    urhs[i][j] = urhs[i][j]
			-v[i][j] * (-u[i][j+2] + 8.0 * ( u[i][j+1] - u[i][j-1] ) + u[i][j-2]) / (12.0 * dy)
			-Math.abs(v[i][j]) * ( u[i][j+2] - 4.0 * u[i][j+1] + 6.0 * u[i][j] - 4.0 * u[i][j-1] + u[i][j-2] ) / ( 4.0 * dy);
		    vrhs[i][j] = vrhs[i][j]
			-v[i][j] * (-v[i][j+2] + 8.0 * ( v[i][j+1] - v[i][j-1] ) + v[i][j-2]) / (12.0 * dy)
			-Math.abs(v[i][j]) * ( v[i][j+2] - 4.0 * v[i][j+1] + 6.0 * v[i][j] - 4.0 * v[i][j-1] + v[i][j-2] ) / ( 4.0 * dy);
		}
	    }
	}		    

	// (5) 更新
	for (int i = 1; i < mx+1; i++) {
	    for (int j = 1; j < my+1; j++) {
	        if (outbox(i,j)) { //物体の内部は計算しない
		    u[i][j] = u[i][j] + dt * urhs[i][j];
		    v[i][j] = v[i][j] + dt * vrhs[i][j];
		}
	    }
	}
    }


	// 箱の外であればTrueを返す
	public static boolean outbox (int i, int j) {
		if ((i1 <= i && i <= i2) && (j1 <= j && j <= j2)) {
			return false;
		} else if ((i3 <= i && i <= i4) && (j3 <= j && j <= j4)) {
			return false;
		} else {
			return true;
		}
	}





    // 実行部
    public static void main (String[] args) throws java.lang.Exception
    {
	long start = System.currentTimeMillis();
	setflw();
	setgrd();
	slvflw();
	long end = System.currentTimeMillis();
	System.out.println("Run Time: " + (end-start) + "[ms]" );
    }
}
