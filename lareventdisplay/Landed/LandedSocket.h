#ifndef LANDEDSOCKET_H
#define LANDEDSOCKET_H

#include "art/Framework/Services/Registry/ServiceMacros.h"

namespace fhicl { class ParameterSet; }
namespace art   { class ActivityRegistry; }
namespace art   { class Worker; }
namespace art   { class InputSource; }
namespace art   { class EventID; }
namespace art   { class Event; }

#include <zlib.h>
#include <boost/asio/local/stream_protocol.hpp>
#include <boost/asio/io_service.hpp>

namespace evd
{

  class LandedSocket
  {
  public:
    LandedSocket(fhicl::ParameterSet const&, art::ActivityRegistry& reg);

  private:
    void postBeginJob();
    void postBeginJobWorkers(art::InputSource*,
			     std::vector<art::Worker*> const&);
    void preProcessEvent(art::Event const&);
    void postProcessEvent(art::Event const&); 

  private:
    art::InputSource* inputSource_;
    boost::asio::local::stream_protocol::endpoint* endpoint_;
    boost::asio::local::stream_protocol::socket* socket_;
    boost::asio::io_service service_;
    unsigned long record_;

  public:
    bool connect(void);
    void sendGeometry(std::string);
    void sendEvent(int, int, int, int, int);
    void sendHit(double, double, double, double, double, int);
    void sendVertex(double, double, double, double, int, int);
  };

}

DECLARE_ART_SERVICE(evd::LandedSocket, LEGACY)

#endif //ifndef LANDEDSOCKET_H
