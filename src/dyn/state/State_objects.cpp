#include "../State_objects.h"

void (State::*init_velocities1)(double)      = &State::init_velocities;
void (State::*init_velocities2)(double,VECTOR,VECTOR) = &State::init_velocities;


void export_State_objects(){

  class_<State>("State",init<>())
      .def(init<System&>())
      .def(init<System&,Thermostat&>())
      .def(init<System&,Barostat&>())
      .def(init<System&,Thermostat&,Barostat&>())

      .def_readwrite("E_kin",&State::E_kin)
      .def_readwrite("E_kin_tr",&State::E_kin_tr)
      .def_readwrite("E_kin_rot",&State::E_kin_rot)
      .def_readwrite("E_pot",&State::E_pot)
      .def_readwrite("E_tot",&State::E_tot)
      .def_readwrite("H",&State::H)
      .def_readwrite("curr_T",&State::curr_T)
      .def_readwrite("curr_P",&State::curr_P)
      .def_readwrite("curr_V",&State::curr_V)
      .def_readwrite("H_NP",&State::H_NP)
      .def_readwrite("L_tot",&State::L_tot)
      .def_readwrite("P_tot",&State::P_tot)


      .def("set",&State::set)
      .def("show_info",&State::show_info)

      .def("update",&State::update)
      .def("cool",&State::cool)
      .def("init_velocities",init_velocities1)
      .def("init_velocities",init_velocities2)

      .def("set_md",&State::set_md)
      .def("init_md",&State::init_md)
      .def("run_md",&State::run_md)
      .def("run_md_respa",&State::run_md_respa)
  ;

}


