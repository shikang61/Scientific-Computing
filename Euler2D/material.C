#include "material.H"

Material::Material(const std::string &type, const double &gamma, const double &p_ref, const bool &MHD)
    : m_type(type), m_gamma(gamma), m_p_ref(p_ref), m_MHD(MHD) {};

std::string Material::type() const
{
    return m_type;
};

double Material::gamma() const
{
    return m_gamma;
};

double Material::p_ref() const
{
    return m_p_ref;
};

bool Material::MHD() const
{
    return m_MHD;
};
