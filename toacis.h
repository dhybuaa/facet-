

#define TRIMEDGE_01	0x0000000000000001		//三角片的第一条边或者四边片的下边界在原始曲面的裁剪边上
#define TRIMEDGE_02	0x0000000000000010		//三角片的第一条边或者四边片的右边界在原始曲面的裁剪边上
#define	TRIMEDGE_03	0x0000000000000100		//三角片的第一条边或者四边片的上边界在原始曲面的裁剪边上
#define TRIMEDGE_04	0x0000000000001000		//三角片的第一条边或者四边片的左边界在原始曲面的裁剪边上
#define LINJIEPATTRI	0x0000000000010000	// 三角片或者四边片是否是临界三角片或临界四边片，这类三角片或者四边片的特点是：
// 角点对应曲面上点的法矢 z 分量都 <= 0




//基本数学常数和精度
#ifndef PI
#define PI		(2*acos(0))
#endif
#ifndef PI2
#define PI2		(4*acos(0))
#endif
#define PAI  3.14159265358979323
#define PAI2 6.28318530717958646

#ifndef ZERO
#define ZERO		1.0e-5
#endif
#ifndef MATHZERO
#define MATHZERO	1.0e-12
#endif
#ifndef KNOTZERO
#define KNOTZERO	1.0e-10
#endif
#ifndef ANGLEZERO
#define ANGLEZERO	1.0e-5
#endif
#ifndef DISTZERO
#define DISTZERO	1.0e-6
#endif
#ifndef AREAZERO
#define AREAZERO	1.0e-4
#endif
#ifndef LINKZERO
#define LINKZERO	2.0e-3
#endif
#ifndef TOACIS_H_INCLUDED
#define TOACIS_H_INCLUDED
struct curve;
//typedef OFFSET   OFFSETT;

typedef  enum {
	ACISPOINT,	//点
	ACISLINE,	//线段
	ACISARC,		//圆弧
	ACISPOLYLINE,//折线
	ACISSPLINE,	//样条线
	ACISOFFSET,	//等距线
	ACISCURCMP,	//组合线，与ME原有的组合线概念不同
	//以下三个曲线类型为兼容ME的函数接口而设置,它们不会存在于系统链中,
	//请不要在外部使用.
	ACISCURVE,	//ME的单曲线
	ACISCOMPOSITE,
	ACISCURVECMP,	//ME的组合线
	ACISNONE,	//超类，无类型
} CURVETYPE;

typedef struct xy {
	double 		x;
	double 		y;
}XY, * PXY, ** PPXY;



typedef struct xyz {
public:
	double 		x;
	double 		y;
	double 		z;

}XYZ, * PXYZ, ** PPXYZ;



typedef struct xyzw {
public:
	double 		x;
	double 		y;
	double 		z;
	double      w;

}XYZW, * PXYZW, ** PPXYZW;



typedef struct box2d {

	double		xmin, ymin;
	double		xmax, ymax;
}BOX2D, * PBOX2D;



typedef struct box3d {

	double		xmin, ymin, zmin;
	double		xmax, ymax, zmax;
}BOX3D, * PBOX3D;

typedef struct pot {
	XYZ p;

}POT, * PPOT;

typedef struct line {
	XYZ p0;
	XYZ p1;
}LINE, * PLINE;


typedef struct polyline {
	int n;
	PXYZ p;
}POLYLINE, * PPOLYLINE;



typedef struct mespline {

	char   nurbsflag; //NURBS,BSPLINE
	char   plnflag;   //NOPLANE,ISPLANE
	double a, b, c, d;
	int    n;
	int    degree;
	PXYZW  pw;
	double* t;
}MESPLINE, * PMESPLINE;


typedef struct arc {
	XYZ     pc, px, py;
	double  r, a;
	PMESPLINE  spline;
}ARC, * PARC;

typedef struct  offset {

	double		dis;
	XYZ		refdir;
	struct  curve* oldcur;
	struct  curve* newcur;
}OFFSET, * POFFSET;

typedef union primitive
{
	POT* point;
	LINE* line;
	ARC* arc;
	POLYLINE* polyline;
	MESPLINE* spline;
	OFFSET* offset;
} PRIMITIVE, * PPRIMITIVE;


typedef struct  mecurve {

	CURVETYPE    type;
	PRIMITIVE    element;
	struct mecurve* next;
}MECURVE, * PMECURVE;



typedef  enum
{
	//SRULE,
	//SROTATE,
	SMESH,
	//SOFFSET,
	//SPATH,
	//SADDCUR,
	//SADDVER,
	//ACISSURF
} SURTYPE;

typedef struct  surrule {

	PMECURVE  curve1;
	PMECURVE  curve2;

}SURRULE, * PSURRULE;


typedef struct surrotate {

	PMECURVE      curve;
	PLINE       axis;
	double      sa;
	double      ea;

}SURROTATE, * PSURROTATE;


typedef struct surmesh {

	int     nurbsflag;         //1----non rational;0----rational.
	int     udegree, vdegree;
	int     unum, vnum;
	PPXYZW  vert;
	double* u, * w;//uv节点矢量数组？
}SURMESH, * PSURMESH;


typedef union surprimitive
{
	//SURRULE* rule;
	//SURROTATE* rotate;
	SURMESH* mesh;
	//class GeomSurOff* offset;
	//class NcSurPath* surpath;		
	//class NcSurAddcur* curpath;		
	//class NcSurAddpt* ptpath;		
} SURPRIMITIVE, * PSURPRIMITIVE;


typedef struct suroffset {

	SURTYPE type;
	double dis;
	SURPRIMITIVE oldsur;
}SUROFFSET, * PSUROFFSET;



typedef struct str2D {
	int    n;
	XY* p;
	char* flag;
	struct str2D* next;
} STR2D, * PSTR2D;
#define SIZESTR2D sizeof(STR2D)

typedef struct str3D {
	int  n;
	XYZ* p;
	struct str3D* next;
} STR3D, * PSTR3D;
#define SIZESTR3D sizeof(STR3D)

typedef struct composite {

	char   breakflag; //NOBREAK,BREAK
	char   plnflag; //NOPLANE,ISPLANE
	double a, b, c, d;
	PMECURVE curve;
	struct  composite* next;
}COMPOSITE, * PCOMPOSITE;

typedef struct curvecmp {

	//unsigned long	handle;
//	int			color;
	//char		layer;
	//char		blank;
	//char		pick_id;
	PMECURVE		curve;
	PCOMPOSITE	cmp;
	PMECURVE		des_cur;
	BOX3D		closebox;

	struct curvecmp* next;
}CURVECMP, * PCURVECMP;


typedef struct mesurf {

	//unsigned long	handle;
	SURTYPE			type;
	//int				shadecolor;
	//int				color;
	//char 			layer;
	//char			blank;
	//char			pick_id;
	char			trimicon; // trimicon==1 be trimed;
	BOX3D			closebox;
	//UWSTR* des_pt;
	SURPRIMITIVE	element;
	//PCURVECMP		spaceoutloop;
	//PCURVECMP		spaceinloop;
	PSTR2D			outloop;
	PSTR2D			inloop;

	//int				uline, wline;
	//PMECURVE			disp_cur;

	struct  mesurf* next;
	//PHANDL			hdlink;

	//PSURTRI			trilink;
}SURF, * PSURF;


//Cut from inter.h by HuangLing.1998.11.24.
typedef struct triangle {
	int  flag;	// 用来标记三角形的边是否是曲面的边界
	// 缺省值是 0x000000000000，用 MallocTriangle() 函数得到的三角片中赋这个初始值
	// 注意：运算的时候使用“宏”，不要使用数值，因为“宏”代表的数值可能随时会发生变化。
	// 操作请用“与(&)”、“或(|)”等。

	// 低 1 - 4 位表示三角片的裁剪信息，对于三角片只用到了 3 位而已
	// 用宏TRIMEDGE_01来表示，0x0000000000000001	表明 p0-p1 是曲面的边界
	// 用宏TRIMEDGE_02来表示，0x0000000000000010	表明 p1-p2 是曲面的边界
	// 用宏TRIMEDGE_03来表示，0x0000000000000100	表明 p2-p0 是曲面的边界

	// 低 5 位表示三角片是否是临界三角片，这类三角片的特点是：
	// 三个参数角点所对应曲面上点的法矢 z 分量全部 <= 0
	// 用宏LINJIEPATTRI来表示，0x0000000000010000
	XYZ  p0;
	XYZ  p1;
	XYZ  p2;
	XY   para0;
	XY   para1;
	XY   para2;
	XYZ  v0;
	XYZ  v1;
	XYZ  v2;
	struct triangle* link;
} *PTRIANGLE, TRIANGLE;
#define SIZETRIANGLE  sizeof(TRIANGLE)

typedef struct patch {
	int  flag;	//用来标记四边片的边是否是曲面的边界
	// 缺省值是 0x000000000000，用 MallocPatch() 函数得到的四边片中赋这个初始值

	// 注意：运算的时候使用“宏”，不要使用数值，因为“宏”代表的数值可能随时会发生变化。
	// 操作请用“与(&)”、“或(|)”等。

	// 低 1 - 4 位表示四边片的裁剪信息
	// 用宏TRIMEDGE_01来表示，0x0000000000000001	表明 p0-p1 是曲面的边界
	// 用宏TRIMEDGE_02来表示，0x0000000000000010	表明 p1-p2 是曲面的边界
	// 用宏TRIMEDGE_03来表示，0x0000000000000100	表明 p2-p3 是曲面的边界
	// 用宏TRIMEDGE_04来表示，0x0000000000001000	表明 p3-p0 是曲面的边界

	// 低 5 位表示四边片是否是临界四边片，这类四边片的特点是：
	// 四个参数角点所对应曲面上点的法矢 z 分量全部 <= 0
	// 用宏LINJIEPATTRI来表示，0x0000000000010000
	XYZ  p0;
	XYZ  p1;
	XYZ  p2;
	XYZ  p3;
	XY   para1;
	XY   para2;
	struct patch* link;
} *PPATCH, PATCH;
#define SIZEPATCH  sizeof(PATCH)


#define DSCPATCH PATCH
#define PDSCPATCH PPATCH
#define DSCTRIANGLE TRIANGLE
#define PDSCTRIANGLE PTRIANGLE


typedef struct linestr
{
	LINE       line;
	short      fg;			// 在求离散边界的时候，这个标记的用法如下：
	// 1：表示这个边界是正常的顺序
	// -1: 表示这个边界反过来了
	XY         uw0, uw1;
	PSURF      sur;
	struct linestr* next;
}LINESTR, * PLINESTR;





#ifndef _INC_GEDSCSUR_H
#include <fstream>
using namespace std;


/*------------------------------------------------------------------
定义宏
------------------------------------------------------------------*/

//参数子环点标记值
//ORITRIMPOINT：  表明该点与链表中下一点形成的直线属于原始曲面的参数定义环边的一部分
//NOTORITRIMPOINT：表明该点与链表中下一点形成的直线不属于原始曲面的参数环边
#define ORITRIMPOINT	-1	//在曲面参数环上的点，也包括因分裂产生的新点，该新点在参数环上。
#define NOTORITRIMPOINT	1   //该点不在参数环上

//交点重复度值
//注意：当交点的重复度值大于NOTCROSS时，交点的方向标志值才起作用
#define BADCROSS	-1
#define NOTCROSS	0
#define SINGLECROSS 1
#define DOUBLECROSS 2

//交点的方向标记值
//注意：在经过裁剪后得到的参数子环中，该标记表明某一交点的下一个交点
//      是应该沿当前的参数父环寻找，还是沿当前的分裂直线寻找
#define TOLOOP		-1
#define TOEDGE		1
#define NODIR		0

//参数子域的类型值
//OUTTYPE：参数外子域
//INTYPE：参数内子域
//TRIMTYPE：参数裁剪子域
#define OUTTYPE		0
#define INTYPE		1
#define TRIMTYPE	2

//参数内子域的属性值
//PUREIN：纯参数内子域，即仅用四个交点描述
//UNPUREIN：非纯参数内子域，即子域的边界上存在补缝点
#define PUREIN		0
#define UNPUREIN	1

//分裂直线与参数子环相交得到的交点的属性值
//TRIMCROSS：若交点为该值，则表明包含该交点的参数子域类型值为TRIMTYPE
//ONCROSS：若交点为该值，则不能立刻判断包含该交点的参数子域类型值是否为TRIMTYPE
#define TRIMCROSS	0
#define ONCROSS		1
#define EDGECROSS	2

//二维点的旋转特性值
//定义：设二维点为cur, 其上一点为pre, 下一点为next;
//      按照pre, cur, next依次排列的顺序形成的旋转方向称为cur的旋转特性值

// 特别注意：以下三个宏的值不能够更改，因为程序中利用了 -1, 0, 1 三个数字的数字特点
#define CLOCKWISE		-1	//顺时针
#define STRAIGHT		0	//共线
#define ANTICLOCKWISE	1	//逆时针

//位置标志
// 特别注意：以下三个宏的值不能够更改，因为程序中利用了 -2, -1, 0, 1, 2 三个数字的数字特点
#define	DSCLEFT			-2	// 左
#define DSCLEFTON		-1	// 左端点
#define	DSCMID			0	// 中间点
#define	DSCRIGHTON		1	// 右端点
#define	DSCRIGHT		2	// 右

//曲面形状标记值
#define ISTWIST	0	//曲面扭曲
#define NOTFLAT 1	//曲面不“平”
#define ISFLAT	2	//曲面满足精度要求

//曲线形状标记值
#define NOTSTRAIGHT	0	//曲线不“直”
#define ISSTRAIGHT	1	//曲线满足精度要求

//曲面分裂方向标记值
#define UDIRECTION	0	//U向分裂
#define VDIRECTION	1	//V向分裂

//曲线的类型标记值
#define DSCPOINT	0	//曲线为点
#define DSCLINE		1	//曲线为直线
#define	DSCARC		2	//曲线为圆弧
#define	DSCSPLINE	3	//曲线为样条或等距线

//离散方式宏
#define PATDSCMODE	0	//四边片混合三角片离散方式
#define TRIDSCMODE	1	//三角片离散方式

//离散要求
#define	CRACKDSC	0	//有缝离散
#define NOCRACKDSC	1	//无缝离散

//环标志
#define NOVALIDLOOP		0	//没有环
#define VALIDLOOP		1	//有环

//交点类型
#define NO_CROSSTYPE	0	// 无属性
#define OUT_CROSSTYPE	1	// 与外环的交点
#define	IN_CROSSTYPE	-1	// 与内环的交点

//环记录方法
#define LOOP_NUM_VAL	0	// 用点的个数来记录环的个数
#define	LOOP_PTER_VAL	1	// 用指针来记录环的个数

//离散裁剪环的优化标记
#define	LOOP_NOOPTIMAL	0	// 裁剪环不进行优化
#define LOOP_OPTIMAL	1	// 裁剪环进行优化

//子域中环的好坏标志
#define LOOP_VALID_OK		0		// 好环
#define LOOP_VALID_FAIL		-1		// 坏环

//内子域(包括纯内子域和非纯内子域)的四条边界的裁剪属性
//即那些边界位于裁剪环上, 以便于得到原始曲面的边界线.
//本系统确定outloop为NULL, 代表该子域为内子域, 
//而inloop的地址值表示该子域的裁剪属性, 如果inloop为NULL, 

/* 以下宏移动到 math3d.h 中

#define TRIMEDGE_01		1	// (0x0000000000000001)		三角片的第一条边或者四边片的下边界在原始曲面的裁剪边上
#define TRIMEDGE_02		3	// (0x0000000000000010)		三角片的第一条边或者四边片的右边界在原始曲面的裁剪边上
#define	TRIMEDGE_03		7	// (0x0000000000000100)		三角片的第一条边或者四边片的上边界在原始曲面的裁剪边上
#define TRIMEDGE_04		15	// (0x0000000000001000)		三角片的第一条边或者四边片的左边界在原始曲面的裁剪边上

*/

//则表明内子域的四条边界都不在原始曲面的边界线上
//该值可以进行"与", "或"及"非"运算. 
//!!!!!注意该前提是outloop为NULL

/*------------------------------------------------------------------
定义宏结束
------------------------------------------------------------------*/

/*------------------------------------------------------------------
定义数据结构
------------------------------------------------------------------*/

// 裁剪节点
// 在离散，分割的过程中，参数域上会不断地出现新的点。
// 每个点，内存中仅有一份。
// used ++ 就是用来记录，外部引用这个点的次数。
typedef struct trimNode
{
	unsigned int id;	// 该节点的标号，在中间计算过程中，该标号还有其他的临时作用
	int used;			// 该节点当前被引用的个数
	XY	p;				// 参数值

	int	 original;// original = 1，表示该点是曲面裁剪环上的原始参数点，否则就不是// bjt-2024

	struct trimNode* pre;
	struct trimNode* next;
}TRIMNODE, * PTRIMNODE, ** PPTRIMNODE;

typedef struct dsctvector		// 用来记录一个线段的长度，以及单位矢量
{
	double len;				// 线段长度
	XY	st;					// 单位矢量
}DSCTVECTOR, * PDSCTVECTOR, ** PPDSCTVECTOR;

//参数子环点
// 每个domain中，如果有裁剪环，那么裁剪环中，记录的就是 trimPoint 的顺序列表
// 而 trimPoint 实际上仅仅是引用包含了 trimNode
typedef struct trimPoint
{
	int			flag; // 该标记用来记录裁剪环上该点，是不是属于最开始曲面裁剪环上的点，如果是就标记为 ORITRIMPOINT
	PDSCTVECTOR	tv;
	PTRIMNODE	ptrim;
	struct trimPoint* pre, * next;
}TRIMPOINT, * PTRIMPOINT, ** PPTRIMPOINT;

//离散拓扑图点
// 该结构，就是为了 domain 的 corner 准备的，实际上，他的内容就是一个 trimNode
// 但需要记住角点的上下左右关系，才有这个结构
typedef struct tolPoint		//(tol:topological)
{
	PTRIMNODE	ptrim;
	struct tolPoint* left, * right, * up, * down;
}TOLPOINT, * PTOLPOINT, ** PPTOLPOINT;

//三维空间子环点
typedef struct hitrimPoint		//(hi:high dimension)
{
	XYZ p;
	struct hitrimPoint* pre, * next;
}HITRIMPOINT, * PHITRIMPOINT, ** PPHITRIMPOINT;

//交点
typedef struct crossPoint
{
	int	type;				//记录是与外环的交点还是与内环的交点，暂时未用
	int leftdir, rightdir;	//交点分别在左、右参数子域中的方向标记(the direction of crosspoint in left or right domain)
	int leftdeg, rightdeg;	//交点分别在左、右参数子域中的重复度(the degeree of crosspoint in left or right domain)
	int fleft, gleft;		//交点在左参数子域中的子环点标记(the flag of crosspoint in left domain)
	int fright, gright;		//交点在右参数子域中的子环点标记(the flag of crosspoint in right domain)
	PTRIMPOINT	new_point;	//指向新生成的裁剪点
	PTRIMPOINT	old_point;		//指向原来的裁剪点
	PTRIMPOINT	nextpoint;	//指向下一个裁剪点
	struct crossPoint* pre, * next;
}CROSSPOINT, * PCROSSPOINT, ** PPCROSSPOINT;

//参数子环
typedef struct trimLoop
{
	int clockwise;			// 顺逆时针
	PTRIMPOINT head;		//参数子环的链表头
	struct trimLoop* next;
}TRIMLOOP, * PTRIMLOOP, ** PPTRIMLOOP;

//参数子域
// outloop, inloop，表示该子域的外环和内环
// 如果该环上有原始裁剪环上的点，该点的 flag 就会是 ORITRIMPOINT，如果该点是新产生的点，就不是这个 flag 标记
// 如果子域已经细分到没有任何裁剪环了，那么 outloop / inloop 就为 null，那么细分三角片的时候，就用 corner 的信息
typedef struct dscDomain
{
	int		loop_valid;
	PPTOLPOINT	corner;			//指向参数子域角点的指针
	PTRIMLOOP	outloop;		//指向参数子域中参数外子环的指针
	PTRIMLOOP	inloop;			//指向参数子域中参数内子环的指针
	struct dscDomain* next;
}DSCDOMAIN, * PDSCDOMAIN, ** PPDSCDOMAIN;

//曲面离散全程信息
typedef struct dscInf	//(Inf:Information)
{
	int	flag;				//标志离散过程是否需要优化曲面的裁剪环
	int maxid;				// 标号指示器，用于赋值 trimNode 的 id，这样每个 trimNode 的 id 都是唯一的
	int foutloop, finloop;	//离散曲面的内外环标志
	double maxdis;			//离散精度
	fstream* fxfs;		//文件指针
	PSURF surf;				//原始曲面指针
	PDSCDOMAIN	domhead;	//曲面离散参数子域链表的“虚”链头(the virtual head of domain link)
	PTOLPOINT	tolhead;	//离散拓扑图的链头(the head of topological point)
	PTRIMNODE	trmhead;	//裁剪节点链表
}DSCINF, * PDSCINF, ** PPDSCINF;

// Double 数链表结构
typedef struct double_link
{
	double para;
	struct double_link* next;
}DUBLINK, * PDUBLINK, ** PPDUBLINK;

class fx_Triangle
{
public:
	int	flag;
	unsigned int p0, p1, p2;
	class fx_Triangle* next;
	fx_Triangle() { flag = 0; p0 = p1 = p2 = 0; next = NULL; }
};

class fx_Patch
{
public:
	int flag;
	unsigned int p0, p1, p2, p3;
	class fx_Patch* next;
	fx_Patch() { flag = 0; p0 = p1 = p2 = p3 = 0; next = NULL; }
};
#define MAX_HASH 10007 
typedef struct {
	bool used;
	XY data;
} HashEntry;
typedef class fx_Triangle	FXTRIANGLE;
typedef class fx_Triangle* PFXTRIANGLE;
typedef class fx_Patch		FXPATCH;
typedef class fx_Patch* PFXPATCH;

/*------------------------------------------------------------------
定义数据结构结束
------------------------------------------------------------------*/
#endif 

#define _INC_GEDSCSUR_H
#endif

/*------------------------------------------------------------------
本文件结束(The End)
------------------------------------------------------------------*/