/*********************************************************************
 * @file   facet.hxx
 * @brief
 * @details
 * @author  dhy
 * @date    2023.6.17
 *********************************************************************/
#pragma once

class outcome;
class ENTITY_LIST;
class AcisOptions;
class BODY;
// 新 demo
outcome aei_FACET_DEMO_DHYTEST(ENTITY_LIST& output_ents, AcisOptions* ptrAcisOpt = nullptr);//测试案例1，星形体，生成stl文件，路径为d盘
outcome aei_FACET_DEMO_DHYTEST2(ENTITY_LIST& output_ents, AcisOptions* ptrAcisOpt);//测试案例2，简单body
outcome aei_FACET_DEMO_DHYTEST3(ENTITY_LIST& output_ents, AcisOptions* ptrAcisOpt);//测试案例3，带内环的洞体
outcome gme_api_facet_entity(BODY* body);//facet离散核心api函数
int gme_api_create_refinement();//控制离散参数api函数
int set_surface_tol(double tol);//控制曲面离散精度api函数，tol缺省值为0.05
