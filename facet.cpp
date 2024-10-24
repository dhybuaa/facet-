/*********************************************************************
 * @file    facet.cpp
 * @brief
 * @details
 * @author  dhy
 * @date    2023.6.17
 *********************************************************************/

#include "facet.hxx"

#include "access.hpp"
#include "acis/acis.hxx"
#include "acis/alltop.hxx"
#include "acis/ckoutcom.hxx"
#include "acis/unitvec.hxx"
#include "acis/logical.h"
#include "acis/position.hxx"
#include "acis/bs3surf.hxx"
#include "acis/sps3srtn.hxx"
#include "acis/spldef.hxx"
#include  "acis/cstrapi.hxx" 
#include "acis/kernapi.hxx" 
#include "acis/lists.hxx" 
#include "acis/vector.hxx"
#include "acis/allsurf.hxx"
#include "acis/allsfdef.hxx"
#include "acis/fileinfo.hxx"
#include "acis/position.hxx"
#include "acis/param.hxx"
#include <acis/sp3srtn.hxx>
#include "acis/bipoly.hxx"
#include "acis/intcurve.hxx"
#include <acis/getbox.hxx>
#include "acis/sgquery.hxx"
#include <acis/sp3crtn.hxx>
#include "acis/bnd_crv.hxx"
#include <acis/sps3crtn.hxx>
#include "toacis.h"
#include "acis/intcurve.hxx"
#include <acis/getbox.hxx>
#include "acis/sgquery.hxx"
#include <acis/sp3crtn.hxx>
#include "acis/bnd_crv.hxx"
#include <acis/sps3crtn.hxx>
#include "template_simple_api.hxx"
#include <stdio.h>
 // 自定义头文件
#include "template_simple_api.hxx"
#include "acis/unitvec.hxx"
// ACIS
#include "acis/logical.h"
#include "acis/acis.hxx"
#include "acis/api.err"
#include "acis/api.hxx"
#include "acis/body.hxx"
#include "acis/check.hxx"
#include "acis/cstrapi.hxx"
#include "acis/module.hxx"
#include "acis/primtive.hxx"
#include "stdafx.h"
#include "toacis.h"
#include <fstream>
#include <stdio.h>
#include "acis/sps3srtn.hxx"
#include "acis/spldef.hxx"
#include  "acis/cstrapi.hxx" 
#include "acis/kernapi.hxx" 
#include "acis/lists.hxx" 
#include "acis/vector.hxx"
#include "acis/allsurf.hxx"
#include "acis/allsfdef.hxx"
#include "acis/fileinfo.hxx"
#include "acis/position.hxx"
#include "acis/param.hxx"
#include <acis/sp3srtn.hxx>
#include "acis/bipoly.hxx"
#include "acis/intcurve.hxx"
#include <acis/getbox.hxx>
#include "acis/sgquery.hxx"
#include <acis/sp3crtn.hxx>
#include "acis/bnd_crv.hxx"
#include <acis/sps3crtn.hxx>
#include "acis/boolapi.hxx"
#include "acis/sfsfint.hxx"
#include "acis/fn2.hxx"
#include "acis/ssi_inf.hxx"
#include "acis/skin_opts.hxx"
#include "acis/skin.hxx"
#include "acis/skin_intr.hxx"
#include "acis/skin_opts.hxx"
#include "acis/skinapi.hxx"
#include "acis/surdef.hxx"
#include "acis/surf_utl.hxx"
#include "acis/surface.hxx"
#include "acis/trace.hxx"
#include "acis/transf.hxx"
#include "acis/wire_utl.hxx"
#include "acis/geolib.hxx"
#include <acis/ptfcenum.hxx>
#include <omp.h>



#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif




int g_PatTriMaxLenControl = 0;	// 离散是否判断四边片或者三角片的最大边长
double g_PatTriMaxLen = 5;		// 定义离散中四边片或者三角片的边长值
int g_Body_dscflag = 0;			// 实体是否已经离散了

double facet_tol;//曲面离散精度


#define ABSOLUTEZERO	1.0e-12		//绝对零(The zero of equation)
#define PARAMETERZERO	1.0e-10		//参数零(The zero for parameters)
#define PARAANGLEZERO	1.0e-10		//角度零(The zero of angle for parameter line)
#define JDGRELATIONZERO	1.0e-6		//判断两条直线关系时用到的零概念
#define DSCSAMEZERO		5.0e-6		//初始裁剪环线段的最小长度(The minimum span for parameter)
#define REMOVEPARAZERO	5.0e-4		//删除初始裁剪点的零值
#define DIMENSIONZERO	1.0e-8		//空间零(The zero for three dimension elements)
#define XYZCORLINEZERO	1.0e-10		//共线控制(The zero for xyz_corline)
#define XYCORLINEZERO	1.0e-10		//共线控制(The zero for xy_corline)
#define MINDIVZERO		5.0e-6		//最小参数域控制(The zero for minimum division)
#define XIACHANGCEOF	3			//空间四边片狭长系数
#define NORMALANGLE		10			//控制法矢之间夹角，单位是“度”

int GeomDscSurf_d2_d3(PSURF surf, double maxdis, int fcrack, int fpatch, PDSCTRIANGLE* tri)
{
	//PTRIFACET tri;
	return 1;
}



/*------------------------------------------------------------------
函数原形声明
-----------------------------------------------------;-------------*/

/*------------------------------------------------------------------
本文件结束(The End)
------------------------------------------------------------------*/

/*
	函数名称：	_GeomDscTestSurfLoop_
	函数描述：	为了测试不同裁剪情况，这里做一个测试函数，用来修改曲面的裁剪状况

	参数说明：
				dscinf	全局裁剪信息
	返回值:
*/
void mathSwapXYZ(PXYZ m_p0, PXYZ m_p1)
{
	XYZ p;
	p = *m_p0;
	*m_p0 = *m_p1;
	*m_p1 = p;
	return;
}


PSTR2D InitSurfOutloop(void)
{
	PSTR2D str = NULL;
	PXY    p = NULL;
	char* s = NULL;

	str = new STR2D;
	str->n = 4;
	str->p = p = new XY[str->n + 1];
	str->flag = s = new char[str->n + 1];

	p[0].x = 0.; p[0].y = 0; s[0] = 0;
	p[1].x = 1.; p[1].y = 0.; s[1] = 0;
	p[2].x = 1.; p[2].y = 1.; s[2] = 0;
	p[3].x = 0.; p[3].y = 1.; s[3] = 0;
	p[4].x = 0.; p[4].y = 0.; s[4] = 0;
	return(str);
}
void _GeomDscTestSurfLoop_(PDSCINF dscinf)
{
	return;

}
void AntiXYZString(PXYZ p, int n)
{
	PXYZ low = NULL;
	int i;


	low = p;
	for (i = 0; i < (n + 1) / 2; i++)
		mathSwapXYZ((low + i), (low + n - i));

	return;
}
void FreeDscTriangleLinkChain(PDSCTRIANGLE tri)
{
	PDSCTRIANGLE next;
	while (NULL != tri)
	{
		next = tri->link;
		delete(tri);
		tri = next;
	}
	return;
}
//DHY7-23
void my_convert_line(PMESPLINE pMesh, bs3_curve bs3)
{
	int dim, deg, num_ctrlpts, num_knots;
	logical rat;
	SPAposition* ctrlpts = NULL;
	double* weights = NULL, * knots = NULL;
	bs3_curve_to_array(bs3, dim, deg, rat, num_ctrlpts, ctrlpts, weights, num_knots, knots);
	pMesh->nurbsflag = 11111;//待填入
	pMesh->plnflag = 11111;//待填入
	pMesh->degree = dim;//阶数？次数？
	pMesh->n = num_ctrlpts;
	pMesh->pw = new XYZW[pMesh->n];
	pMesh->t = knots;
	for (int i = 0; i < num_ctrlpts - 1; ++i) {
		pMesh->pw[i].x = ctrlpts[i].coordinate(0);
		pMesh->pw[i].y = ctrlpts[i].coordinate(1);
		pMesh->pw[i].z = ctrlpts[i].coordinate(2);
	}
}
//DHY7-23
void my_convert_box(SPApar_box myBoxPtr, bs3_surface bs3, PBOX2D ACISabox)
{

	myBoxPtr = bs3_surface_range(bs3);
	SPAinterval ub = bs3_surface_range_u(bs3);
	SPAinterval vb = bs3_surface_range_v(bs3);

	vb.start_pt();
	vb.end_pt();
	ub.start_pt();
	ub.end_pt();


	int i = bs3_surface_nspans_u(bs3);
	int j = bs3_surface_nspans_v(bs3);
	//bs3_surface_bispan_poly(0, 0, bs3);
	//SPAinterval ub1 = bs3_surface_range_u(bs3);
	//SPAinterval vb1 = bs3_surface_range_v(bs3);


}

//董浩宇修改，加入全局变量metocurve
bs3_curve metocurve;
int GeneratorSpline(PMESPLINE curve, double t, PXYZ o_p)
{


	SPAinterval rangcur = bs3_curve_range(metocurve);//获取定义域的取值范围
	double u = rangcur.start_pt() + t * (rangcur.end_pt() - rangcur.start_pt());
	SPAposition acis_point = bs3_curve_position(u, metocurve);
	o_p->x = acis_point.coordinate(0);//x
	o_p->y = acis_point.coordinate(1);//y
	o_p->z = acis_point.coordinate(2);//z
	return 1;
}


//acis曲线转换me曲线7.26g
//晚上修改，未对其做归一化处理
PMESPLINE SplineAcisToMe(bs3_curve aciscurve)
{

	PMESPLINE curve = new MESPLINE;
	int dim, deg, num_ctrlpts, num_knots;
	logical rat;
	SPAposition* ctrlpts = ACIS_NEW SPAposition();
	double* weights = ACIS_NEW double;
	double* knots = ACIS_NEW double;
	bs3_curve_to_array(aciscurve, dim, deg, rat, num_ctrlpts, ctrlpts, weights, num_knots, knots);

	if (weights == NULL)
		curve->nurbsflag = 0;//待填入
	else
		curve->nurbsflag = 1;


	curve->plnflag = 0;//待填入
	curve->degree = deg;//
	curve->n = num_ctrlpts - 1;
	curve->pw = new XYZW[curve->n + 1];
	//curve->t = knots;
	curve->t = new double[curve->n + curve->degree + 2];
	SPAinterval range = bs3_curve_range(aciscurve);

	for (int i = 0; i < num_knots; i++) {

		//curve->t[i] = (curve->t[i] - range.start_pt()) / (range.end_pt() - range.start_pt());
		curve->t[i] = (knots[i] - range.start_pt()) / (range.end_pt() - range.start_pt());
	}


	for (int i = 0; i < num_ctrlpts; ++i) {
		curve->pw[i].x = ctrlpts[i].coordinate(0);
		curve->pw[i].y = ctrlpts[i].coordinate(1);
		curve->pw[i].z = ctrlpts[i].coordinate(2);
		if (curve->nurbsflag == 0)
			curve->pw[i].w = 1.;
		else
			curve->pw[i].w = weights[i];
	}

	metocurve = aciscurve;

	ACIS_DELETE STD_CAST knots;
	ACIS_DELETE STD_CAST weights;
	ACIS_DELETE ctrlpts;
	return curve;
}
//下面是闫光荣做的按误差离散曲线的函数
int  DescrtSplineByTol(PMESPLINE i_spline, double maxdis, PXYZ* o_p)
{
	int divnum = 0;



	return divnum;
}



void SaveTrianglesToAscSTL(char* filename, PDSCTRIANGLE i_tri)
{

	PDSCTRIANGLE tri;
	FILE* fp = NULL;

	if (NULL == i_tri) return;

	fp = fopen(filename, "w+t");

	fprintf(fp, "solid ascii\n");
	for (tri = i_tri; tri != NULL; tri = tri->link)
	{

		fprintf(fp, "  facet normal %.6Le %.6Le %.6Le\n", tri->v0.x, tri->v0.y, tri->v0.z);
		//	fprintf(fp,"  facet normal %.5lf %.5lf %.5lf\n",tri->vec0.x,tri->vec0.y,tri->vec0.z);
		fprintf(fp, "    outer loop\n");
		//fprintf(fp, "      vertex %.5lf %.5lf %.5lf\n", tri->p0.x, tri->p0.y, tri->p0.z); 2022-12-22
		fprintf(fp, "      vertex %.6Le %.6Le %.6Le\n", tri->p0.x, tri->p0.y, tri->p0.z);
		//fprintf(fp, "      vertex %.5lf %.5lf %.5lf\n", tri->p1.x, tri->p1.y, tri->p1.z);2022-12-22
		fprintf(fp, "      vertex %.6Le %.6Le %.6Le\n", tri->p1.x, tri->p1.y, tri->p1.z);
		//fprintf(fp, "      vertex %.5lf %.5lf %.5lf\n", tri->p2.x, tri->p2.y, tri->p2.z);
		fprintf(fp, "      vertex %.6Le %.6Le %.6Le\n", tri->p2.x, tri->p2.y, tri->p2.z);
		fprintf(fp, "    endloop\n");
		fprintf(fp, "  endfacet\n");
	}
	fprintf(fp, "endsolid\n");
	fclose(fp);
	return;
}

void mapToUnitSquare(SPApar_box myBoxPtr, double& x, double& y)
{
	SPAinterval ub = myBoxPtr.u_range();
	SPAinterval vb = myBoxPtr.v_range();
	x = (x - ub.start_pt()) / (ub.end_pt() - ub.start_pt());
	y = (y - vb.start_pt()) / (vb.end_pt() - vb.start_pt());
}
double mapValue(double& x, double& u, double& v) {
	return u + (v - u) * x;
}


//acis曲面转me曲面
PSURMESH mesh_acis_to_me(bs3_surface bs3)
{
	PSURMESH pMesh = new SURMESH;
	double* weights = ACIS_NEW double;//权因子数组不是二维数组！

	SPAposition* ctrlpts = ACIS_NEW SPAposition();//控制顶点数组
	int nurbsflag = 1;


	if (bs3_surface_rational_u(bs3) == 0 && bs3_surface_rational_v(bs3) == 0)	//ture 1，false 0
		pMesh->nurbsflag = 1;
	else pMesh->nurbsflag = 0;
	pMesh->udegree = bs3_surface_degree_u(bs3);//阶数？次数？
	pMesh->vdegree = bs3_surface_degree_v(bs3);
	//pMesh->unum = bs3_surface_ncu(bs3) - 1;
	//pMesh->vnum = bs3_surface_ncv(bs3) - 1;
	bs3_surface_weights(bs3, pMesh->unum, pMesh->vnum, weights);
	if (weights == NULL)
	{
		nurbsflag = 0;
	}
	bs3_surface_control_points(bs3, pMesh->unum, pMesh->vnum, ctrlpts);
	pMesh->unum -= 1;
	pMesh->vnum -= 1;
	int 	num_knots_u;
	int 	num_knots_v;
	double* uknots = ACIS_NEW double;
	double* vknots = ACIS_NEW double;
	bs3_surface_knots_u(bs3, num_knots_u, uknots);
	bs3_surface_knots_v(bs3, num_knots_v, vknots);
	pMesh->u = new double[num_knots_u];
	pMesh->w = new double[num_knots_v];

	SPApar_box myBoxPtr;

	myBoxPtr = bs3_surface_range(bs3);
	SPAinterval ub = myBoxPtr.u_range();
	SPAinterval vb = myBoxPtr.v_range();
	for (int i = 0; i < num_knots_u; i++) {
		pMesh->u[i] = uknots[i];
		pMesh->u[i] = (pMesh->u[i] - ub.start_pt()) / (ub.end_pt() - ub.start_pt());
	}
	for (int i = 0; i < num_knots_v; i++) {
		pMesh->w[i] = vknots[i];
		pMesh->w[i] = (pMesh->w[i] - vb.start_pt()) / (vb.end_pt() - vb.start_pt());
	}

	int ctu = pMesh->unum;
	int ctv = pMesh->vnum;
	pMesh->vert = new PXYZW[ctu + 1];
	for (int i = 0; i <= ctu; i++) {
		pMesh->vert[i] = new XYZW[ctv + 1];
	}
	for (int i = 0; i <= ctu; i++) {
		for (int j = 0; j <= ctv; j++) {
			int index = i * (ctv + 1) + j;
			pMesh->vert[i][j].x = ctrlpts[index].coordinate(0);//x
			pMesh->vert[i][j].y = ctrlpts[index].coordinate(1);//y
			pMesh->vert[i][j].z = ctrlpts[index].coordinate(2);//z

			//为保证测试案例2通过做修改7.23
			// 
			if (nurbsflag == 1)
				pMesh->vert[i][j].w = weights[index];//w	
			else
				pMesh->vert[i][j].w = 1.;//w
		}
	}
	ACIS_DELETE STD_CAST uknots;
	ACIS_DELETE STD_CAST vknots;
	ACIS_DELETE STD_CAST weights;
	ACIS_DELETE ctrlpts;
	return pMesh;//
}
//按照给定数据构建曲线
void my_construct_bs3curve(bs3_curve& bs3)
{
	int degree = 4;
	int rational = FALSE;         //是否为有理曲线?否
	int  closed = FALSE;             //是否为闭合曲线?否
	int  periodic = FALSE;             //是否为周期曲线?否
	int num_ctrlpts = 6;
	SPAposition ctrlpts[6] = { SPAposition(40, 0, -40), SPAposition(40, 5, -20), SPAposition(40, 10,-40),
		SPAposition(40, 15,-20), SPAposition(40, 20,-40), SPAposition(40, 25, -20) };
	double weights[6] = { 1.0, 1.0, 1.0,
						  1.0, 1.0, 1.0 };
	double ctrlpt_tol = SPAresabs;
	int   num_knots = 11;//节点向量中节点的个数
	double   knots[11] = { 0,0,0,0,0,0.5,1,1,1,1,1 };//节点向量
	double knot_tol = SPAresabs;

	bs3 = bs3_curve_from_ctrlpts(
		degree,
		rational,
		closed,
		periodic,
		num_ctrlpts,
		ctrlpts,
		weights,
		ctrlpt_tol,
		num_knots,
		knots,
		knot_tol);
	PMESPLINE mecur = SplineAcisToMe(bs3);
	PXYZ despoint = NULL;
	int n = DescrtSplineByTol(mecur, 0.01, &despoint);

}
//按照给定数据构建曲面
void my_construct_bs3(bs3_surface& bs3)
{
	// Hard code values for this example.
	// Construct a wiggly quarter pipe.
	int degree_u = 3;
	int degree_v = 3;
	logical rational_u = FALSE;
	logical rational_v = TRUE;
	int form_u = 0;  // open
	int form_v = 0;  // open
	int pole_u = 0;  // no singularities at u bounds
	int pole_v = 0;  // no singularities at v bounds
	int num_ctrlpts_u = 4;
	int num_ctrlpts_v = 4;
	// Values are stored in a 1-D array,
	// with v values varying first.
	SPAposition ctrlpts[16];
	// column 0
	/*
	*	ctrlpts[0] = SPAposition(0.0, 0.0, 0.0);
	ctrlpts[1] = SPAposition(10.0, 5.0, 10.0);
	ctrlpts[2] = SPAposition(20.0, 5.0, 10.0);
	ctrlpts[3] = SPAposition(30.0, 0.0, 0.0);
	// column 1
	ctrlpts[4] = SPAposition(0.0, 10.0, 0.0);
	ctrlpts[5] = SPAposition(10.0, 15.0, 10.0);
	ctrlpts[6] = SPAposition(20.0, 15.0, 10.0);
	ctrlpts[7] = SPAposition(30.0, 10.0, 0.0);
	// column 2
	ctrlpts[8] = SPAposition(0.0, 20.0, 0.0);
	ctrlpts[9] = SPAposition(10.0, 25.0, 10.0);
	ctrlpts[10] = SPAposition(20.0, 25.0, 10.0);
	ctrlpts[11] = SPAposition(30.0, 20.0, 0.0);
	// column 3
	ctrlpts[12] = SPAposition(0.0, 30.0, 0.0);
	ctrlpts[13] = SPAposition(10.0, 35.0, 10.0);
	ctrlpts[14] = SPAposition(20.0, 35.0, 10.0);
	ctrlpts[15] = SPAposition(30.0, 30.0, 0.0);
	*/
	ctrlpts[0] = SPAposition(0.0, 0.0, 5.0);
	ctrlpts[1] = SPAposition(10.0, 5.0, 10.0);
	ctrlpts[2] = SPAposition(20.0, 5.0, 15.0);
	ctrlpts[3] = SPAposition(30.0, 0.0, 1.0);
	// column 1
	ctrlpts[4] = SPAposition(1.0, 10.0, 6.0);
	ctrlpts[5] = SPAposition(14.0, 15.0, 13.0);
	ctrlpts[6] = SPAposition(20.0, 15.0, 10.0);
	ctrlpts[7] = SPAposition(32.0, 10.0, 2.0);
	// column 2
	ctrlpts[8] = SPAposition(3.0, 20.0, 8.0);
	ctrlpts[9] = SPAposition(12.0, 25.0, 10.0);
	ctrlpts[10] = SPAposition(21.0, 25.0, 14.0);
	ctrlpts[11] = SPAposition(34.0, 20.0, 3.0);
	// column 3
	ctrlpts[12] = SPAposition(4.0, 30.0, 10.0);
	ctrlpts[13] = SPAposition(13.0, 35.0, 13.0);
	ctrlpts[14] = SPAposition(25.0, 35.0, 17.0);
	ctrlpts[15] = SPAposition(36.0, 32.0, 0.0);
	// column 4

	// Values are stored in a 1-D array,
	// with same order as ctrlpts.
	double weights[16] = { 1.0, 1, 1.0,
						  1.0, 1, 1.0,
						  1.0, 1, 1.0,
						  1.0, 1, 1.0,
						  1.0, 1, 1.0,1 };
	double ctrlpt_tol = SPAresabs;

#if 1
	// Use end knot multiplicity of degree.
	int num_knots_u = 8;
	double knots_u[8] = { 0.0, 0.0, 0.0, 0,3,3,3,3 };
	int num_knots_v = 8;
	double knots_v[8] = { 1.0, 1.0, 1.0, 1,4, 4,4,4 };
#else
	// Use end knot multiplicity equal to the (degree+1).
	int num_knots_u = 9;
	double knots_u[9] = { 0.0, 0.0, 0.0, 0.0, 2.0, 4.0, 4.0, 4.0, 4.0 };
	int num_knots_v = 6;
	double knots_v[6] = { 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 };
#endif
	double knot_tol = SPAresabs;

	// Construct the b-spline surface.
	bs3 = bs3_surface_from_ctrlpts(
		degree_u, rational_u, form_u, pole_u, num_ctrlpts_u,
		degree_v, rational_v, form_v, pole_v, num_ctrlpts_v,
		ctrlpts, weights, ctrlpt_tol,
		num_knots_u, knots_u, num_knots_v, knots_v, knot_tol);
	printf("ddddd");
}

//bs3转face
void my_construct_face(bs3_surface bs3, FACE*& f)
{
	// Construct a spline surface from the bs3_surface.
	spline* spl = ACIS_NEW spline(bs3);

	// Construct a face from the spline surface.
	// This will construct its own copy of the spline surface.
	api_make_face_from_surface(spl, f);

	// Delete our local copy of the spline surface.
	if (spl != NULL)
		ACIS_DELETE spl;
}
//获取曲面公共边7/24dhy
void getinsect(FACE*& bsbool1, FACE*& bsbool2, help_point*& list)
{
	//FACE* f = ACIS_NEW FACE();
	//my_construct_face(bsbool1, f);

//	help_point* list;
	const SPAbox& box = get_face_box(bsbool1);
	const surface& sf1 = bsbool1->geometry()->equation();
	const surface& sf2 = bsbool2->geometry()->equation();
	surf_surf_int* boolen = d3_sf_sf_int(
		sf1,
		sf2,
		box,
		SPAresfit,
		list,
		ALL_CURVES,
		FALSE,
		0.0,
		0.0,
		FALSE
	);
}
//在acis中生成stl文件操作
int  save_acis_stl(PSURF sur)
{

	PDSCTRIANGLE tri = NULL;
	GeomDscSurf_d2_d3(sur, facet_tol, 1, 1, &tri);
	char filemename[256];
	strcpy(filemename, "D:\\test.stl");
	SaveTrianglesToAscSTL(filemename, tri);
	return 1;

}
//dhy创建，EDGE转换成bs3curve
bs3_curve Acis_Edge_to_bs3(EDGE*& my_edge)
{

	CURVE* my_curve = my_edge->geometry();
	const  curve& dhy_curve = my_curve->equation();
	curve* dhy_curve2 = dhy_curve.copy_curve();
	bs3_curve bs3curve1 = NULL;
	if (dhy_curve2->type() == intcurve_type) {
		// Retrieving the bs3_curve only if the curve type is 'intcurve'
		bs3curve1 = ((intcurve*)dhy_curve2)->cur();


	}
	//int ctpoint = bs3_curve_num_ctlpts(bs3curve1);

	//出现问题，在demo1调用时未正确输出，输出edgelist【2】为空7/30,进行修改
	SPAparameter start = my_edge->start_param();
	SPAparameter end = my_edge->end_param();
	double& actual_tol = SpaAcis::NullObj::get_double();
	bs3_curve bs3curve2 = bs3_curve_make_cur(dhy_curve, start, end, 0, *(double*)NULL_REF);

	if (bs3curve1 != NULL)
		return bs3curve1;
	else
		return bs3curve2;
}


//dhy创建，FACE转换成bs3surface
bs3_surface Acis_FACE_to_bs3(FACE*& my_face)
{

	SURFACE* fs = my_face->geometry();
	const surface& sfs = fs->equation();
	SPAbox re = get_face_box(my_face);
	SPAtransf& tr = sg_get_transform(fs);
	double& actual_fit = SpaAcis::NullObj::get_double();
	SPApar_transf& pt = SpaAcis::NullObj::get_par_transf();
	bs3_surface bs3faceacis = bs3_surface_make_sur(sfs, re, 0, *(double*)NULL_REF, *(SPApar_transf*)NULL_REF);
	return bs3faceacis;
}
PSURF MeshToSurface(PSURMESH mesh)
{
	PSURF sur = NULL;

	sur = new SURF;
	sur->type = SMESH;
	sur->element.mesh = mesh;
	sur->outloop = InitSurfOutloop();
	sur->inloop = NULL;
	sur->trimicon = 0;

	return sur;
}

#define MAX_HASH_SIZE 10007
bool hash_set[MAX_HASH_SIZE] = { false }; // 初始化哈希集合为全部false

int hash_function(XY elem) {
	long long hash = (long long)(elem.x * 1000) % MAX_HASH_SIZE;
	hash = hash * 1000 + (long long)(elem.y * 1000) % MAX_HASH_SIZE;
	return (int)hash % MAX_HASH_SIZE;
}

void remove_duplicates(XY arr[], int* size) {
	if (*size <= 1) return;

	int write_index = 0;

	for (int read_index = 0; read_index < *size; read_index++) {
		XY elem = arr[read_index];
		int hash = hash_function(elem);

		if (!hash_set[hash]) {
			hash_set[hash] = true;
			arr[write_index] = elem;
			write_index++;
		}
	}

	*size = write_index;
	memset(hash_set, false, sizeof(hash_set));
}
//为acisfacetome添加的结构体
typedef struct parapt {
	XY p;
	char flag;
	struct parapt* next;
}Parapt, * PParapt;
double mathGetPtPtDis(PXYZ p0, PXYZ p1)
{
	double dx, dy, dz;
	dx = p1->x - p0->x;
	dy = p1->y - p0->y;
	dz = p1->z - p0->z;
	return sqrt(dx * dx + dy * dy + dz * dz);
}
bs3_surface metoacis_sur;// 定义一个全局变量，用于存储acis的surface
//重要函数,outloop与inloop
PSURF AcisFACEtoMeSurfinloop(FACE* f, double tol)
{
	bs3_surface bs3 = Acis_FACE_to_bs3(f);
	PSURMESH pmesh = mesh_acis_to_me(bs3);
	//metoacis_sur = bs3_surface_copy(bs3dhy);//全局变量赋值
	metoacis_sur = bs3;
	PSURF sur = MeshToSurface(pmesh);
	PParapt hstr = NULL, str = NULL;



	//边生成
	int i = 0;
	ENTITY_LIST edge_list, in_edge_list;
	SPAposition s0, e0, s1, e1;
	int num_lpt = 0;
	int n;//= 800; // 初始数组大小
	//char* flag = (char*)malloc(n * sizeof(char));
	//XY* strpt = (XY*)malloc(n * sizeof(XY));
	char* flag = NULL;
	XY* strpt = NULL;
	PXYZ o_p = NULL;

	//outloop操作
	bs3_curve testcurve;
	int anti = 0;
	COEDGE* coedge = NULL;
	XYZ startp, endp;
	int outloopflag = 0;
	int outloopcount = 0;
	sur->inloop = NULL;
	loop_type type;
	int 	info[2];



	for (LOOP* loop = f->loop(); loop; loop = loop->next())
	{
		api_loop_type(loop, type, info);
		int j = loop->size();
		for (coedge = loop->start(), i = 0; i < loop->size(); coedge = coedge->next(), i++) {
			if (type == loop_periphery) {              //如果是outloop
#pragma omp critical
				edge_list.add(coedge->edge());
			}                                         //有的有向边里面没有edge，此处i循环了112次得到四几何边
		}
	}
#pragma omp parallel for schedule(dynamic)

	char filemename[256];

	int num_edges = edge_list.count();
	/*for (int i = 0; i < num_edges; i++) {

		EDGE* my_edge = (EDGE*)edge_list[i];
		testcurve = Acis_Edge_to_bs3(my_edge);
		PMESPLINE spline = SplineAcisToMe(testcurve);
		{
			strcpy(filemename, "D:\\2.txt");
			FILE* fp = fopen(filemename, "w+t");
			for (int k = 0; k < spline->n; k++)
			{
				fprintf(fp, "%f\n", spline->pw[k].z);
			}
			fclose(fp);
		}
		if (spline->pw->z == 6)
			printf("6");
		int n = DescrtSplineByTol(spline, tol, &o_p);
		SPAinterval ub = bs3_surface_range_u(bs3);
		SPAinterval vb = bs3_surface_range_v(bs3);


		if (i != 0)
		{
			anti = 0;

			if (mathGetPtPtDis(&endp, &o_p[0]) > 0.01) anti = 1;
			if (anti == 1)
			{
				AntiXYZString(o_p, n);
			}

		}

		startp = o_p[0];
		endp = o_p[n];
		for (int j = 0; j < n; j++)
		{
			SPAposition acispt(o_p[j].x, o_p[j].y, o_p[j].z);
			SPApar_pos acispos = bs3_surface_invert(acispt, bs3);//给定x，y，z得到u，v
			//strpt[num_lpt].x = (acispos.u - ub.start_pt()) / (ub.end_pt() - ub.start_pt());
			//strpt[num_lpt].y = (acispos.v - vb.start_pt()) / (vb.end_pt() - vb.start_pt());
			if (hstr == NULL)
				hstr = str = (PParapt)malloc(sizeof(Parapt));
			else
				str = str->next = (PParapt)malloc(sizeof(Parapt));
			str->next = NULL;
			str->p.x = (acispos.u - ub.start_pt()) / (ub.end_pt() - ub.start_pt());
			str->p.y = (acispos.v - vb.start_pt()) / (vb.end_pt() - vb.start_pt());
			if (j == 0)
				//flag[num_lpt] = 0;
				str->flag = 0;
			else
				//flag[num_lpt] = 108;
				str->flag = 108;
			num_lpt++;
		}
	}
	flag = (char*)malloc((num_lpt + 1) * sizeof(char)); //该参数环的参数点个数为num_lpt+1。
	strpt = (XY*)malloc((num_lpt + 1) * sizeof(XY));
	num_lpt = 0;
	for (str = hstr; str != NULL; str = str->next)
	{
		strpt[num_lpt].x = str->p.x;
		strpt[num_lpt].y = str->p.y;
		flag[num_lpt] = str->flag;
		num_lpt++;
	}
	strpt[num_lpt].x = strpt[0].x;
	strpt[num_lpt].y = strpt[0].y;
	flag[num_lpt] = 0;

	while (hstr != NULL)
	{
		str = hstr->next;
		free(hstr);
		hstr = str;
	}
	//此处应判断环的方向，确保外环为逆时针，内环为顺时针。
	sur->outloop = new STR2D;
	sur->trimicon = 1;
	sur->outloop->n = num_lpt;
	sur->outloop->p = strpt;
	sur->outloop->flag = flag;
	sur->outloop->next = 0;


	//inloop
	for (LOOP* loop = f->loop(); loop; loop = loop->next())
	{
		//每次更新inloop dhy

		hstr = NULL;
		str = NULL;
		api_loop_type(loop, type, info);
		int innum_lpt = 0;
		char* inflag = NULL;//(char*)malloc(n * sizeof(char));
		XY* instrpt = NULL;//(XY*)malloc(n * sizeof(XY));
		PXYZ ino_p = NULL;
		int j = loop->size();
		in_edge_list.clear();

		for (coedge = loop->start(), i = 0; i < loop->size(); coedge = coedge->next(), i++) {
			if (type == loop_hole)           //如果是inloop
				in_edge_list.add(coedge->edge());

		}
#pragma omp parallel for schedule(dynamic)
		int in_num_edges = in_edge_list.count();
		for (int i = 0; i < in_num_edges; i++) {
			EDGE* my_edge = (EDGE*)in_edge_list[i];
			testcurve = Acis_Edge_to_bs3(my_edge);
			PMESPLINE spline = SplineAcisToMe(testcurve);
			int n = DescrtSplineByTol(spline, tol, &ino_p);
			SPAinterval ub = bs3_surface_range_u(bs3);
			SPAinterval vb = bs3_surface_range_v(bs3);
			if (i != 0)
			{
				anti = 0;

				if (mathGetPtPtDis(&endp, &ino_p[0]) > 0.01) anti = 1;
				if (anti == 1)
				{
					AntiXYZString(ino_p, n);
				}

			}
			startp = ino_p[0];
			endp = ino_p[n];
			for (int j = 0; j < n; j++)
			{
				SPAposition acispt(ino_p[j].x, ino_p[j].y, ino_p[j].z);
				SPApar_pos acispos = bs3_surface_invert(acispt, bs3);//给定x，y，z得到u，v
				//instrpt[innum_lpt].x = (acispos.u - ub.start_pt()) / (ub.end_pt() - ub.start_pt());
				//instrpt[innum_lpt].y = (acispos.v - vb.start_pt()) / (vb.end_pt() - vb.start_pt());
				if (hstr == NULL)
					hstr = str = (PParapt)malloc(sizeof(Parapt));
				else
					str = str->next = (PParapt)malloc(sizeof(Parapt));
				str->next = NULL;
				str->p.x = (acispos.u - ub.start_pt()) / (ub.end_pt() - ub.start_pt());
				str->p.y = (acispos.v - vb.start_pt()) / (vb.end_pt() - vb.start_pt());
				if (j == 0)//innum_lpt == 0)
					str->flag = 0;
				else
					str->flag = 108;
				innum_lpt++;
			}
		}
		//	instrpt[innum_lpt].x = instrpt[0].x;
		//	instrpt[innum_lpt].y = instrpt[0].y;
		//	inflag[innum_lpt] = 0;
			//inloop链表操作
		inflag = (char*)malloc((innum_lpt + 1) * sizeof(char)); //该参数环的参数点个数为num_lpt+1。
		instrpt = (XY*)malloc((innum_lpt + 1) * sizeof(XY));
		innum_lpt = 0;
		for (str = hstr; str != NULL; str = str->next)
		{
			instrpt[innum_lpt].x = str->p.x;
			instrpt[innum_lpt].y = str->p.y;
			inflag[innum_lpt] = str->flag;
			innum_lpt++;
		}
		instrpt[innum_lpt].x = instrpt[0].x;
		instrpt[innum_lpt].y = instrpt[0].y;
		inflag[innum_lpt] = 0;

		while (hstr != NULL)
		{
			str = hstr->next;
			free(hstr);
			hstr = str;
		}

		STR2D* new_inloop = new STR2D;
		new_inloop->n = innum_lpt;
		new_inloop->p = instrpt;
		new_inloop->flag = inflag;
		new_inloop->next = sur->inloop;
		sur->inloop = new_inloop;

	}

	{
		char filemename[256];
		strcpy(filemename, "D:\\1.txt");
		FILE* fp = fopen(filemename, "w+t");
		for (int k = 0; k < num_lpt; k++)
		{
			//fprintf(fp, "%f %f\n", strpt[k].x, strpt[k].y);

			fprintf(fp, "%f\n", strpt[k].x);
		}
		fclose(fp);
	}
	{
		char filemename[256];
		strcpy(filemename, "D:\\2.txt");
		FILE* fp = fopen(filemename, "w+t");
		for (int k = 0; k < num_lpt; k++)
		{
			fprintf(fp, "%f\n", strpt[k].y);
		}
		fclose(fp);
	}

	*/

	return sur;
}

//两曲面拼接

outcome gme_api_facet_entity(BODY* body) {
	/* Initialization Block */
	BODY* skin_curved_star_body = NULL;
	API_BEGIN
		PDSCTRIANGLE htri = NULL, tri = NULL, stri = NULL;
	for (LUMP* lump = body->lump(); lump; lump = lump->next()) {
		// 遍历 BODY 中的 LUMP
		// 处理当前 LUMP
		for (SHELL* shell = lump->shell(); shell; shell = shell->next()) {
			// 遍历 LUMP中的 SHELL
			// 处理当前 SHELL
			for (FACE* face = shell->face(); face; face = face->next()) {

				{
					PSURF sur = AcisFACEtoMeSurfinloop(face, facet_tol);
					stri = NULL;
					//sur->trimicon = 1;
					GeomDscSurf_d2_d3(sur, facet_tol, 1, 1, &stri);
					if (htri == NULL)htri = stri;
					else
					{
						tri = htri;
						while (tri->link != NULL)tri = tri->link;
						tri->link = stri;
					}
					//	FreeSurf_in_ACIS(sur);	//break;
				}

			}
		}
	}
	char filemename[256];
	strcpy(filemename, "D:\\test.stl");
	SaveTrianglesToAscSTL(filemename, htri);
	API_END
		return result;
}

outcome aei_FACET_DEMO_DHYTEST2(ENTITY_LIST& output_ents, AcisOptions* ptrAcisOpt) {
	/* Initialization Block */
	FACE* f = NULL;
	FACE* contf = NULL;
	/* API Call Block */
	API_BEGIN
		bs3_surface bs3;
	SPApar_box myBoxPtr;
	my_construct_bs3(bs3);
	my_construct_face(bs3, f);
	gme_api_create_refinement();

	int i = 0;
	ENTITY_LIST edge_list;
	COEDGE* coedge = NULL;
	for (LOOP* loop = f->loop(); loop; loop = loop->next())
	{
		for (coedge = loop->start(), i = 0; i < loop->size(); coedge = coedge->next(), i++) {

			edge_list.add(coedge->edge());
			//下面部分挪动到此处。

			//SPAposition start=edge->start()->geometry();
		}
	}



	bs3_curve contcurve;
	bs3_curve bs3curvecont;
	EDGE* my_edge = (EDGE*)edge_list[2];
	bs3curvecont = Acis_Edge_to_bs3(my_edge);
	my_construct_bs3curve(contcurve);
	bs3_surface contbs3 = bs3_surface_ruled(contcurve, bs3curvecont);
	my_construct_face(contbs3, contf);


	ENTITY_LIST face_list;
	face_list.add(contf);
	face_list.add(f);
	int num_faces = face_list.count();
	PDSCTRIANGLE htri = NULL, tri = NULL, stri = NULL;
	for (int i = 0; i < num_faces; i++) {
		{
			FACE* intf = (FACE*)face_list[i];
			PSURF sur = AcisFACEtoMeSurfinloop(intf, facet_tol);
			stri = NULL;
			//sur->trimicon = 1;
			GeomDscSurf_d2_d3(sur, facet_tol, 1, 1, &stri);
			if (htri == NULL)htri = stri;
			else
			{
				tri = htri;
				while (tri->link != NULL)tri = tri->link;
				tri->link = stri;
			}
		}

		//break;

	}

	char filemename[256];
	strcpy(filemename, "D:\\test.stl");
	SaveTrianglesToAscSTL(filemename, htri);
	API_END

		if (result.ok()) {
			//ENTITY* ent = nullptr;
			//acis_api_restore_entity("C:\\Users\\yangu\\Desktop\\cuboid.sat", ent);
			output_ents.add(f);
			output_ents.add(contf);
		}
	return result;
}
//简单body，设计内环,并行计算
outcome aei_FACET_DEMO_DHYTEST3(ENTITY_LIST& output_ents, AcisOptions* ptrAcisOpt) {
	/* Initialization Block */
	BODY* block, * pointer;
	//BODY* body = NULL;  // Pointer to tool body
	//FACE* f = NULL;
	/* API Call Block */
	BODY* cyl2;
	BODY* cylinder = NULL;
	API_BEGIN
		api_start_modeller(0);
	api_initialize_kernel();
	BODY* body[2];
	ENTITY_LIST body_list;
	double	width = 50, depth = 50, height = 10;
	gme_api_create_refinement();

	int          low_v = 1;//样条曲面的v 向为s 形状
	int         high_v = 1;
	int         low_u = 2;//样条曲面的u 向为凸起的形状
	int         high_u = 2;
	api_wiggle(width, depth, height, low_v, high_v, low_u, high_u, TRUE, block);
	BODY* cyl = NULL;
	check_outcome(api_make_frustum(60.0, 5.0, 5.0, 5.0, cyl));
	check_outcome(api_subtract(cyl, block));//bool运算

	SPAposition pos1(13.0, 13.0, 8.0);
	api_solid_sphere(pos1,
		6,
		cyl2
	);
	check_outcome(api_subtract(cyl2, block));//bool运算
	SHELL* eg = block->lump()->shell();
	FACE* ff;
	PDSCTRIANGLE htri = NULL, tri = NULL, stri = NULL;
	int facenum = 0;
	SPAposition bottom(0, 0, 0);
	SPAposition top(0.5, 0, 0);
	double maj_rad = 1.0;
	double min_rad = 1.0;



	api_solid_cylinder_cone(bottom, top, maj_rad, min_rad, maj_rad, NULL, cylinder);//一个圆柱

	for (LUMP* lump = block->lump(); lump; lump = lump->next()) {
		// 遍历 BODY 中的 LUMP
		// 处理当前 LUMP
#pragma omp parallel for schedule(dynamic)
		for (SHELL* shell = lump->shell(); shell; shell = shell->next()) {
			// 遍历 LUMP中的 SHELL
			// 处理当前 SHELL
#pragma omp parallel for schedule(dynamic)
			for (FACE* face = shell->face(); face; face = face->next()) {

				//if (facenum == 2)//问题模型法facenum==0
				{
					PSURF sur = AcisFACEtoMeSurfinloop(face, facet_tol);
					stri = NULL;
					//sur->trimicon = 1;
					GeomDscSurf_d2_d3(sur, facet_tol, 1, 1, &stri);
#pragma omp critical
					{
						if (htri == NULL)htri = stri;
						else
						{
							tri = htri;
							while (tri->link != NULL)tri = tri->link;
							tri->link = stri;
						}
					}//break;
				}
				facenum++;
			}
		}
	}
	char filemename[256];
	strcpy(filemename, "D:\\test.stl");
	SaveTrianglesToAscSTL(filemename, htri);
	FreeDscTriangleLinkChain(htri);

	//	acis_api_save_entity("C:\\Users\\yangu\\Desktop\\cuboid.stl", block);
	api_terminate_constructors();
	api_stop_modeller();
	//block = NULL;
	API_END
		output_ents.add(block);

	return result;
}
//复杂body
SPAtransf rot_about_z1(double angle) {
	// SPAvector z_axis(0.0,0.0,1.0);
	return rotate_transf(angle, z_axis);
}
SPAtransf move_along_z1(double distance) {
	SPAvector direction(0.0, 0.0, distance);
	return translate_transf(direction);
}
outcome aei_FACET_DEMO_DHYTEST(ENTITY_LIST& output_ents, AcisOptions* ptrAcisOpt) {
	/* Initialization Block */
	BODY* skin_curved_star_body = NULL;

	//BODY* body = NULL;  // Pointer to tool body
	//FACE* f = NULL;
	/* API Call Block */
	API_BEGIN
		BODY* triangles[12];
	gme_api_create_refinement();




	// Make triangle[0]
	double L = cos(M_PI / 12.0) / sin(M_PI / 12.0);
	SPAposition points[3];
	points[0].set_x(L);
	points[0].set_y(1.0);
	points[0].set_z(0.0);
	points[1].set_x(L + 2.0);
	points[1].set_y(0);
	points[1].set_z(0);
	points[2].set_x(L);
	points[2].set_y(-1);
	points[2].set_z(0.0);
	check_outcome(api_make_wire(NULL, 3, points, triangles[0]));

	// Create triangle[1],...,triangle[11] by copying and rotating triangle[0]
	for (int i = 1; i < 12; i++) {
		check_outcome(api_copy_body(triangles[0], triangles[i]));
		check_outcome(api_transform_entity((ENTITY*)triangles[i], rot_about_z1(i * M_PI / 6.0)));

		// Simultaneously, we also unite the new triangle with triangle[0]
		// s.t. when we are done triangle[0] will be a star.
		check_outcome(api_unite(triangles[i], triangles[0]));
	}

	// We obtain a star
	BODY* star1 = triangles[0];

	// circle1: wire-body
	const SPAposition center(0.0, 0.0, 0.0);
	EDGE* circle1_edge[1];
	check_outcome(api_curve_arc(center, 2, 0.0, 2 * M_PI, circle1_edge[0]));
	check_outcome(api_transform_entity((ENTITY*)circle1_edge[0], move_along_z1(6.0)));

	// Make it into a wire-body
	BODY* circle1 = NULL;
	check_outcome(api_make_ewire(1, circle1_edge, circle1));

	// guide1: EDGE
	EDGE* guide1 = NULL;
	SPAposition guide1_points[3];

	guide1_points[0].set_x(L + 2.0);  // on star1
	guide1_points[0].set_y(0.0);
	guide1_points[0].set_z(0.0);

	guide1_points[1].set_x(6.0);
	guide1_points[1].set_y(0.0);
	guide1_points[1].set_z(3.0);

	guide1_points[2].set_x(2.0);  // on circle1
	guide1_points[2].set_y(0.0);
	guide1_points[2].set_z(6.0);

	check_outcome(api_curve_spline(3, guide1_points, NULL, NULL, guide1));

	// Make an array of wire-bodies
	BODY* wires[2];
	wires[0] = star1;
	wires[1] = circle1;

	// Make an array of guides
	EDGE* guides[1];
	guides[0] = guide1;

	// Set virtualGuides to TRUE
	skin_options opts;
	opts.set_virtualGuides(TRUE);

	// Skin star1 and circle1 with one guide and virtualGuides set to TRUE
	check_outcome(api_skin_wires(2, wires, 1, guides, skin_curved_star_body, &opts, ptrAcisOpt));

	api_del_entity(wires[0]);
	api_del_entity(wires[1]);
	api_del_entity(guide1);
	gme_api_facet_entity(skin_curved_star_body);
	API_END

		if (result.ok()) {
			ENTITY* ent = nullptr;

			output_ents.add(skin_curved_star_body);
		}
	return result;
}
int gme_api_create_refinement() {
	facet_tol = 0.5;
	return 0;
}
int set_surface_tol(double tol) {
	facet_tol = tol;
	return 0;
}