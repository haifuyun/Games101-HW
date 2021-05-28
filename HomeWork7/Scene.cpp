//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    //ergodic scene objects
    for (uint32_t k = 0; k < objects.size(); ++k) {
        //object emit light
        if (objects[k]->hasEmit()){
            //sumbit all emit object area
            emit_area_sum += objects[k]->getArea();
        }
    }
    //random number * all object area
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            //random area < some one object area
            if (p <= emit_area_sum){
                //light sample
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    float pdf_light;
    Intersection lightInter;     
    sampleLight(lightInter,pdf_light);
    
    // TO DO Implement Path Tracing Algorithm here
    //ray intersction
    Intersection inter = intersect(ray);
    

    if (!inter.happened)
    {
        return {0,0,0};
    }
    

    if (inter.m->hasEmission())
    {
        return inter.m->getEmission();
    }


    
    Vector3f blockedCheckDir = lightInter.coords - inter.coords;
    Intersection blockCheckInter = intersect(Ray(inter.coords,blockedCheckDir.normalized()));
    
    Vector3f L_dir = 0.0f;

    if ((blockCheckInter.coords - lightInter.coords).norm() < 0.01)
    {
        Vector3f L_i =  lightInter.emit;
        Vector3f f_r = inter.m->eval(ray.direction,blockedCheckDir.normalized(),inter.normal);
        float cos_theta = dotProduct(blockedCheckDir.normalized(), inter.normal);
        float cos_theta_Light = dotProduct(-blockedCheckDir.normalized(),lightInter.normal);
        float point_light_distance = dotProduct(blockedCheckDir,blockedCheckDir);
        

        L_dir = L_i * f_r * cos_theta * cos_theta_Light / point_light_distance / pdf_light;
    }

    Vector3f L_indir = 0.0f;
    
    //((float)rand() / (RAND_MAX))
    float randomNumber = get_random_float(); 

    if (randomNumber < RussianRoulette)
    {
        Vector3f wi = inter.m->sample(ray.direction,inter.normal).normalized();
        
        Vector3f preL_dir = castRay(Ray(inter.coords,wi),depth);
        Vector3f indir_f_r = inter.m->eval(ray.direction, wi, inter.normal);
        float indir_cos_theta = dotProduct(wi, inter.normal);
        float indir_pdf = inter.m->pdf(ray.direction,  wi, inter.normal);
        L_indir = preL_dir * indir_f_r * indir_cos_theta / indir_pdf / RussianRoulette; 
                
    }
    
    return L_dir + L_indir;
    

}