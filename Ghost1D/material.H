#ifndef MATERIAL_H
#define MATERIAL_H

#include "EoS.H"
#include <sstream>

class Material
{
private:
    std::string m_type;
    double m_gamma;
    double m_p_ref;
    bool m_MHD;

public:
    Material(const std::string &, const double &, const double &, const bool &); // initialisation constructon
    virtual ~Material() = default;

    std::string type() const;
    double gamma() const;
    double p_ref() const;
    bool MHD() const;
};
#endif