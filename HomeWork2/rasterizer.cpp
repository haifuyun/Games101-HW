// clang-format off
//
// Created by goksu on 4/6/19.
//

#include <algorithm>
#include <vector>
#include "rasterizer.hpp"
#include <opencv2/opencv.hpp>
#include <math.h>


rst::pos_buf_id rst::rasterizer::load_positions(const std::vector<Eigen::Vector3f> &positions)
{
    auto id = get_next_id();
    pos_buf.emplace(id, positions);

    return {id};
}

rst::ind_buf_id rst::rasterizer::load_indices(const std::vector<Eigen::Vector3i> &indices)
{
    auto id = get_next_id();
    ind_buf.emplace(id, indices);

    return {id};
}

rst::col_buf_id rst::rasterizer::load_colors(const std::vector<Eigen::Vector3f> &cols)
{
    auto id = get_next_id();
    col_buf.emplace(id, cols);

    return {id};
}

auto to_vec4(const Eigen::Vector3f& v3, float w = 1.0f)
{
    return Vector4f(v3.x(), v3.y(), v3.z(), w);
}

static bool isSameSide(const Eigen::Vector3f& pt1, const Eigen::Vector3f& pt2, 
const Eigen::Vector3f& a, const Eigen::Vector3f& b)
{
// reference: https://blackpawn.com/texts/pointinpoly/default.html
    Eigen::Vector3f cp1 = (b-a).cross(pt1-a);
    Eigen::Vector3f cp2 = (b-a).cross(pt2-a);
if (cp1.dot(cp2)>0)
    {
return true;
    }
else
    {
return false;
    }    
}

static bool insideTriangle(int x, int y, const Vector3f* _v)
{   
    // TODO : Implement this function to check if the point (x, y) is inside the triangle represented by _v[0], _v[1], _v[2]
    // std::clog << "Triangle Point1 " << " " << _v[0] << "\n";
    // std::clog << "Triangle Point2 " << " " << _v[1] << "\n";
    // std::clog << "Triangle Point3 " << " " << _v[2] << "\n";

    Eigen::Vector3f pt;
    pt << (float)x, (float)y, 0.0f;
    bool isInside = isSameSide(pt, _v[0], _v[1], _v[2]) && isSameSide(pt, _v[1], _v[0], _v[2]) && isSameSide(pt, _v[2], _v[0], _v[1]);
    return isInside;

    // int px,py,px1,py1,px2,py2;
    
    // Vector3f point(x,y,1.0f);

    // if (_v[0].x() < _v[1].x())
    // {
    //     if (_v[1].x() < _v[2].x())
    //     {

    //          px = _v[0].x();
    //          py = _v[0].y();
    //          px1 = _v[1].x();
    //          py1 = _v[1].y();
    //          px2 = _v[2].x();
    //          py2 = _v[2].y();

    //     }else{

    //         if (_v[0].x() < _v[2].x())
    //         {
    //             px = _v[0].x();
    //             py = _v[0].y();
    //             px1 = _v[2].x();
    //             py1 = _v[2].y();
    //             px2 = _v[1].x();
    //             py2 = _v[1].y();

    //         }else{

    //             px = _v[2].x();
    //             py = _v[2].y();
    //             px1 = _v[0].x();
    //             py1 = _v[0].y();
    //             px2 = _v[1].x();
    //             py2 = _v[1].y();

    //         }
    //     }
    // }else
    // {
    //     if (_v[1].x() < _v[2].x())
    //     {

    //         if (_v[0].x() < _v[2].x())
    //         {

    //             px = _v[1].x();
    //             py = _v[1].y();
    //             px1 = _v[0].x();
    //             py1 = _v[0].y();
    //             px2 = _v[2].x();
    //             py2 = _v[2].y();

    //         }else{

    //             px = _v[1].x();
    //             py = _v[1].y();
    //             px1 = _v[2].x();
    //             py1 = _v[2].y();
    //             px2 = _v[0].x();
    //             py2 = _v[0].y();

    //         }

    //     }else{

    //         px = _v[2].x();
    //         py = _v[2].y();
    //         px1 = _v[1].x();
    //         py1 = _v[1].y();
    //         px2 = _v[0].x();
    //         py2 = _v[0].y();
            
    //     }
    // }

    // Eigen::Vector3f p1;
    // p1 << px,py,1;
    // Eigen::Vector3f p2;
    // p2 << px1,py2,1;
    // Eigen::Vector3f p3;
    // p3 << px2,py2,1;
    // std::clog << " point 1 \n" << p1 << std::endl;
    // std::clog << " point 2 \n" << p2 << std::endl;
    // std::clog << " point 3 \n" << p3 << std::endl;

    // Eigen::Vector3f v1;
    // v1 = p3 - p1;
    // Eigen::Vector3f v2;
    // v2 = p2 - p3;
    // Eigen::Vector3f v3;
    // v3 = p1 - p2;

    // Eigen::Vector3f v1p;
    // v1p = point - p1;
    // Eigen::Vector3f v2p;
    // v2p = point - p3;
    // Eigen::Vector3f v3p;
    // v3p = point - p2;
    
    // bool inside = v1.cross(v1p).z() > 0 && v2.cross(v2p).z() > 0 && v3.cross(v3p).z() > 0;

    // // std::clog << " v1 c  " << v1.cross(v1p).z() << std::endl;
    // // std::clog << " v2 c  " << v2.cross(v2p).z() << std::endl;
    // // std::clog << " v3 c " << v3.cross(v3p).z() << std::endl;

    // std::clog << " inside " << inside << std::endl;

    // return inside;

}

static std::tuple<float, float, float> computeBarycentric2D(float x, float y, const Vector3f* v)
{
    float c1 = (x*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*y + v[1].x()*v[2].y() - v[2].x()*v[1].y()) / (v[0].x()*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*v[0].y() + v[1].x()*v[2].y() - v[2].x()*v[1].y());
    float c2 = (x*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*y + v[2].x()*v[0].y() - v[0].x()*v[2].y()) / (v[1].x()*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*v[1].y() + v[2].x()*v[0].y() - v[0].x()*v[2].y());
    float c3 = (x*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*y + v[0].x()*v[1].y() - v[1].x()*v[0].y()) / (v[2].x()*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*v[2].y() + v[0].x()*v[1].y() - v[1].x()*v[0].y());
    return {c1,c2,c3};
}

void rst::rasterizer::draw(pos_buf_id pos_buffer, ind_buf_id ind_buffer, col_buf_id col_buffer, Primitive type)
{
    auto& buf = pos_buf[pos_buffer.pos_id];
    auto& ind = ind_buf[ind_buffer.ind_id];
    auto& col = col_buf[col_buffer.col_id];

    float f1 = (50 - 0.1) / 2.0;
    float f2 = (50 + 0.1) / 2.0;

    Eigen::Matrix4f mvp = projection * view * model;
    for (auto& i : ind)
    {
        Triangle t;
        Eigen::Vector4f v[] = {
                mvp * to_vec4(buf[i[0]], 1.0f),
                mvp * to_vec4(buf[i[1]], 1.0f),
                mvp * to_vec4(buf[i[2]], 1.0f)
        };
        //Homogeneous division
        for (auto& vec : v) {
            vec /= vec.w();
        }
        //Viewport transformation
        for (auto & vert : v)
        {
            vert.x() = 0.5*width*(vert.x()+1.0);
            vert.y() = 0.5*height*(vert.y()+1.0);
            vert.z() = -vert.z() * f1 + f2;
        }

        for (int i = 0; i < 3; ++i)
        {
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
        }

        auto col_x = col[i[0]];
        auto col_y = col[i[1]];
        auto col_z = col[i[2]];

        t.setColor(0, col_x[0], col_x[1], col_x[2]);
        t.setColor(1, col_y[0], col_y[1], col_y[2]);
        t.setColor(2, col_z[0], col_z[1], col_z[2]);

        rasterize_triangle(t);
    }
}

//Screen space rasterization
void rst::rasterizer::rasterize_triangle(const Triangle& t) {
    auto v = t.toVector4();
    
    // TODO : Find out the bounding box of current triangle.
    // iterate through the pixel and find if the current pixel is inside the triangle

    // If so, use the following code to get the interpolated z value.
    //auto[alpha, beta, gamma] = computeBarycentric2D(x, y, t.v);
    //float w_reciprocal = 1.0/(alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
    //float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
    //z_interpolated *= w_reciprocal;

    // TODO : set the current pixel (use the set_pixel function) to the color of the triangle (use getColor function) if it should be painted.

    float xMinBox = 0;
    float xMaxBox = 0;
    float yMinBox = 0;
    float yMaxBox = 0;


    
    // auto ext_pt_x = std::minmax_element(v.begin(), v.end(),
    //                                    [](const Eigen::Vector4f& lhs, const Eigen::Vector4f & rhs){
    //                                         return lhs.x() < rhs.x();
    //                                    });
    // auto ext_pt_y = std::minmax_element(v.begin(), v.end(),
    //                                    [](const Eigen::Vector4f& lhs, const Eigen::Vector4f & rhs){
    //                                         return lhs.y() < rhs.y();
    //                                    });
    // xMinBox = floor(ext_pt_x.first->x());
    // xMaxBox = floor(ext_pt_x.second->x());
    // yMinBox = floor(ext_pt_y.first->y());
    // yMaxBox = floor(ext_pt_y.second->y());
    
    //init point array
    int vlength;
    vlength = sizeof(v)/sizeof(Eigen::Vector4f);
    std::clog << "v size " << vlength << '\n';

    Eigen::Vector3f v3f[vlength];
    for (int i = 0; i < vlength; i++)
    {
        std::clog << "v point" << i << " x " << v[i].x() << " y " << v[i].y() << '\n';
        v3f[i] << v[i].x(), v[i].y(), 1.0f;
    }


    xMinBox = v[0].x();
    yMinBox = v[0].y();

    //find min max x y
    for (auto& vec : v)
    {
        if (vec.x() < xMinBox)
        {
            xMinBox = vec.x();
        }
        
        if (vec.x() > xMaxBox)
        {
            xMaxBox = vec.x();
        }
        
        if (vec.y() < yMinBox)
        {
            yMinBox = vec.y();
        }

        if (vec.y() > yMaxBox)
        {
            yMaxBox = vec.y();
        }
    }
    
    /*
    make box range cover triangle
    top high down low
    */
    xMinBox = floor(xMinBox);
    xMaxBox = ceil(xMaxBox);
    yMinBox = floor(yMinBox);
    yMaxBox = ceil(yMaxBox);

    std::clog << xMinBox << " xMinBox" << std::endl;
    std::clog << xMaxBox << " xMaxBox" << std::endl;
    std::clog << yMinBox << " yMinBox" << std::endl;
    std::clog << yMaxBox << " yMaxBox" << std::endl;

    for (int x = xMinBox; x <= xMaxBox; ++x)
    {
        for (int y = yMinBox; y <= yMaxBox; ++y)
        {
            auto[alpha, beta, gamma] = computeBarycentric2D((float)x, (float)y, t.v);
            float w_reciprocal = 1.0/(alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
            float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
            z_interpolated *= w_reciprocal;

            int index = x + width*y;
            
            if(insideTriangle(x,y,v3f) && (z_interpolated < depth_buf[index])){
                Eigen::Vector3f point(x,y,1.0f);
                set_pixel(point,t.getColor());                
                
                depth_buf[index] = z_interpolated;
            };

            // float x_pos = (float)x + 0.5;
            // float y_pos = (float)y + 0.5;
            // auto[alpha, beta, gamma] = computeBarycentric2D(x_pos, y_pos, t.v);
            // float w_reciprocal = 1.0/(alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
            // float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
            // z_interpolated *= w_reciprocal;
            // int buff_ind = x + y*width;
            // if (insideTriangle(x_pos, y_pos, v3f) && (-z_interpolated<depth_buf[buff_ind]))
            // {
            //     set_pixel(Eigen::Vector3f(x, y, 1.0f), t.getColor());
            //     depth_buf[buff_ind] = -z_interpolated;
            // }
        }
    }
}

void rst::rasterizer::set_model(const Eigen::Matrix4f& m)
{
    model = m;
}

void rst::rasterizer::set_view(const Eigen::Matrix4f& v)
{
    view = v;
}

void rst::rasterizer::set_projection(const Eigen::Matrix4f& p)
{
    projection = p;
}

void rst::rasterizer::clear(rst::Buffers buff)
{
    if ((buff & rst::Buffers::Color) == rst::Buffers::Color)
    {
        std::fill(frame_buf.begin(), frame_buf.end(), Eigen::Vector3f{0, 0, 0});
    }
    if ((buff & rst::Buffers::Depth) == rst::Buffers::Depth)
    {
        std::fill(depth_buf.begin(), depth_buf.end(), std::numeric_limits<float>::infinity());
    }
}

rst::rasterizer::rasterizer(int w, int h) : width(w), height(h)
{
    frame_buf.resize(w * h);
    depth_buf.resize(w * h);
}

int rst::rasterizer::get_index(int x, int y)
{
    return (height-1-y)*width + x;
}

void rst::rasterizer::set_pixel(const Eigen::Vector3f& point, const Eigen::Vector3f& color)
{
    //old index: auto ind = point.y() + point.x() * width;
    auto ind = (height-1-point.y())*width + point.x();
    frame_buf[ind] = color;

}

// clang-format on