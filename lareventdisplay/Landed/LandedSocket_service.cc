#include "LandedSocket.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Core/InputSource.h"
#include "art/Framework/Principal/Worker.h"
#include "art/Persistency/Provenance/EventID.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/IO/Root/RootInput.h"

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string.h>

#include <boost/array.hpp>

namespace evd
{
  LandedSocket::LandedSocket(fhicl::ParameterSet const& pset,
		 art::ActivityRegistry& reg):
    endpoint_(nullptr),
    socket_(nullptr)
  {
   
    reg.sPostBeginJob.watch       (this, &LandedSocket::postBeginJob);
    reg.sPostBeginJobWorkers.watch(this, &LandedSocket::postBeginJobWorkers);
    reg.sPreProcessEvent.watch    (this, &LandedSocket::preProcessEvent);
    reg.sPostProcessEvent.watch   (this, &LandedSocket::postProcessEvent);
    
  }

  void LandedSocket::
  postBeginJob()
  {

  }
  void LandedSocket::
  postBeginJobWorkers(art::InputSource* inputsource,
			   std::vector<art::Worker*> const& workers)
  {
    inputSource_=inputsource;
    record_=0;
  }
  void LandedSocket::
  preProcessEvent(art::Event const&)
  {

  }
  void LandedSocket::
  postProcessEvent(art::Event const& event)
  {
    if (socket_!=nullptr)
      {
	art::RootInput* rootinput = dynamic_cast<art::RootInput*>(inputSource_);
	if (!rootinput)
	  {
	    throw cet::exception("Landed") << "input is not root file, therefore not seekable and therefore unsuitable for EventViewer\n";
	  }
	boost::array<char, 128> instr;
	size_t len=socket_->receive_from(boost::asoi::buffer(instr), *endpoint_);
	if (len>0)
	  {
	    long newrec=atol(instr.data());
	    long diff=newrec-record_;
	    diff--;
	    //	    std::cout << "seek: " << diff << std::endl;
	    if (rootinput->seekToEvent(diff, true))
	      {
		std::cout << "seek succeeded\n";
		send(socket_, "OK\n", 3, 0);
		record_=newrec;
	      }
	    else
	      {
		//		std::cout << "no such event, reloading\n";
		if (!rootinput->seekToEvent(event.id(), true))
		  throw cet::exception("Landed") << "problem performing seek on input, lost track of input event\n";
		send(socket_, "NO\n", 3, 0);
		//		std::cout << "reload performed\n";
	      }
	    //wait for message confirmation
	    char conf;
	    if (recv(socket_, &conf, 1, 0)!=1)
	      throw cet::exception("Landed") << "LANDED app did not confirm instruction result\n";	
	    
	  }
	else
	  throw cet::exception("Landed") << "bad instruction from LANDED app\n";
      }
  }

  bool LandedSocket::
  connect(void)
  {
    //connect to landed app
    socket_=socket(AF_UNIX, SOCK_STREAM, 0);
    if (socket_<0)
      {
	std::cerr << "Landed: failed to create socket for connection to LANDED app\n";	
	return false;
      }
    std::ostringstream sockfile;
    sockfile << getenv("HOME") << "/.landed.sock";
    boost::asio::local::stream_protocol::endpoint ep(sockfile.str());
    struct sockaddr_un address;
    memset(&address, 0, sizeof(struct sockaddr_un));
    address.sun_family=AF_UNIX;
    strncpy(address.sun_path, sockfile.str().c_str(), sizeof(address.sun_path));
    address.sun_path[sizeof(address.sun_path)-1]=0;
    if (::connect(socket_, (struct sockaddr*)&address,strlen(address.sun_path)+sizeof(address.sun_family))!=0)
      {
	std::cerr << "Landed: failed to connect to LANDED app\n";
	socket_=-1;
      }
    return socket>=0;
  }

  void LandedSocket::
  sendGeometry(std::string vrml)
  {
    std::ostringstream ss;
    ss << "GEO:"<< vrml.length() << "\n";
    if (send(socket_, ss.str().c_str(), ss.str().length(), 0)!=(ssize_t)ss.str().length() ||
	send(socket_, vrml.c_str(), vrml.length(), 0)!=(ssize_t)(vrml.length()))
      {
	throw cet::exception("Landed") << "failed to send detector geometry to LANDED app\n";	
      }
    char conf;
    if (recv(socket_, &conf, 1, 0)!=1)
      throw cet::exception("Landed") << "LANDED app did not confirm detector geometry\n";
    else
      std::cout << "LANDED app has confirmed detector geometry\n";

  }
  void LandedSocket::
  sendEvent(int nhits, int nvertex, int run, int subrun, int event)
  {
    std::ostringstream ss;
    ss << "EVT:" << record_ << "," << nhits << "," << nvertex << ", " << run << "," << subrun << "," << event << "\n";
    //    std::cout << ss.str();
     if (send(socket_, ss.str().c_str(), ss.str().length(), 0)!=(ssize_t)ss.str().length())
	throw cet::exception("Landed") << "failed to send event introduction to LANDED app\n";	
    char conf;
    if (recv(socket_, &conf, 1, 0)!=1)
      throw cet::exception("Landed") << "LANDED app did not confirm event\n";	
  }
  void LandedSocket::
  sendHit(double x, double y, double z, double e, double ne, int track)
  {
    std::ostringstream ss;
    ss << "HIT:" << x << "," << y << "," << z << "," << e << "," << ne << "," << track << "\n";
    //    std::cout << ss.str();
     if (send(socket_, ss.str().c_str(), ss.str().length(), 0)!=(ssize_t)ss.str().length())
	throw cet::exception("Landed") << "failed to send hit to LANDED app\n";	
    char conf;
    if (recv(socket_, &conf, 1, 0)!=1)
      throw cet::exception("Landed") << "LANDED app did not confirm hit\n";	

  }
  void LandedSocket::
  sendVertex(double x, double y, double z, double e, int id, int pdgcode)
  {
    std::ostringstream ss;
    ss << "VTX:" << x << "," << y << "," << z << "," << e << "," << id << "," << pdgcode << "\n";
     if (send(socket_, ss.str().c_str(), ss.str().length(), 0)!=(ssize_t)ss.str().length())
	throw cet::exception("Landed") << "failed to send vertex to LANDED app\n";	
    char conf;
    if (recv(socket_, &conf, 1, 0)!=1)
      throw cet::exception("Landed") << "LANDED app did not confirm vertex\n";	
  }
}

DEFINE_ART_SERVICE(evd::LandedSocket)
