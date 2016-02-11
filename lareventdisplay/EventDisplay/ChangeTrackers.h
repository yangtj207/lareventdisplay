/**
 * @file   ChangeTrackers.h
 * @brief  Classes detecting configuration changes
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   September 23rd, 2015
 */

#ifndef UTIL_CHANGETRACKERS_H
#define UTIL_CHANGETRACKERS_H

// LArSoft libraries
#include "larcore/SimpleTypesAndConstants/geo_types.h" // geo::PlaneID

// framework libraries
#include "art/Utilities/InputTag.h"
#include "art/Persistency/Provenance/EventID.h"
#include "art/Framework/Principal/Event.h"

// C/C++ standard libraries
#include <string> // std::to_string()
#include <type_traits> // std::enable_if_t


namespace util {
  
  /** **************************************************************************
   * @brief Detects the presence of a new event
   *
   * The state of this class describes the current event by its ID.
   */
  class EventChangeTracker_t {
      public:
    /// Default constructor: no current event, next event is a new one
    EventChangeTracker_t() = default;
    
    /// Constructor: current event as specified
    EventChangeTracker_t(art::Event const& evt): state{evt.id()} {}
    
    /// Constructor: current event as specified by the event ID
    EventChangeTracker_t(art::EventID const& evt_id): state{evt_id} {}
    
    /// @name State query
    /// @{
    /// Returns whether this tracker is in the same state as another
    bool same(EventChangeTracker_t const& as) const
      { return as.eventID() == eventID(); }
    
    /// Returns whether there is a current event
    bool isValid() const { return eventID() != art::EventID(); }
    
    /// Returns whether this tracker is in the same state as another
    bool operator== (EventChangeTracker_t const& as) const { return same(as); }
    
    /// Returns whether this tracker is in a different state than another
    bool operator!= (EventChangeTracker_t const& than) const
      { return !same(than); }
    /// @}
    
    /// @name State change
    /// @{
    /// Forgets the current event
    void clear() { set(art::EventID()); }
    
    /// Sets the current event ID
    void set(art::EventID const& evt_id) { state.event_id = evt_id; }
    
    /// Sets the current event ID from the given event
    void set(art::Event const& evt) { set(evt.id()); }
    
    /// Sets the current event, and returns true if it is changed
    bool update(EventChangeTracker_t const& trk)
      {
        if (same(trk)) return false;
        *this = trk;
        return true;
      }
    /// @}
    
    /// Returns a string representing the current state
    operator std::string() const
      {
        return "R:" + std::to_string(eventID().run())
          + " S:" + std::to_string(eventID().subRun())
          + " E:" + std::to_string(eventID().event());
      }
    
    
      protected:
    /// Returns the current event ID (it might be made public...)
    art::EventID const& eventID() const { return state.event_id; }
    
    
      private:
    struct LocalState_t {
      art::EventID event_id; ///< ID of the current event
    };
    
    LocalState_t state; ///< local state of the tracker (may inherit some more)
    
  }; // EventChangeTracker_t
  
  
  template <typename Stream>
  decltype(auto) operator<<
    (Stream&& out, EventChangeTracker_t const& trk)
    { out << std::string(trk); return std::forward<Stream>(out); }
  
  
  
  /** **************************************************************************
   * @brief Detects the presence of a new event or data product
   *
   * The state of this class describes the current data product input as the
   * event it belongs to (by its ID) and the input tag.
   */
  class DataProductChangeTracker_t: private EventChangeTracker_t {
      public:
    
    /// Default constructor: no current data product
    DataProductChangeTracker_t() = default;
    
    /// Constructor: specifies current event and data product label
    DataProductChangeTracker_t
      (art::Event const& evt, art::InputTag const& label):
      EventChangeTracker_t(evt), state{label}
      {}
    
    /// Constructor: specifies current event ID and data product label
    DataProductChangeTracker_t
      (art::EventID const& evt_id, art::InputTag const& label):
      EventChangeTracker_t(evt_id), state{label}
      {}
    
    
    /// @name State query
    /// @{
    /// Returns the current input label
    art::InputTag const& inputLabel() const { return state.input_label; }
    
    /// Returns whether we are in the same event (the rest could differ)
    bool sameEvent(DataProductChangeTracker_t const& as) const
      { return EventChangeTracker_t::same(as); }
    
    /// Returns whether we have same data product as in "as"
    bool same(DataProductChangeTracker_t const& as) const
      { return sameEvent(as) && (inputLabel() == as.inputLabel()); }
    
    /// Returns whether there is a current event and data product
    bool isValid() const
      {
        return
          EventChangeTracker_t::isValid() && !inputLabel().label().empty();
      }
    
    /// Returns whether the event and input label are the same as in "as"
    bool operator== (DataProductChangeTracker_t const& as) const
      { return same(as); }
    
    /// Returns whether the event or input label are different than in "than"
    bool operator!= (DataProductChangeTracker_t const& than) const
      { return !same(than); }
    /// @}
    
    
    /// @name State change
    /// @{
    /// Forget the current data product
    void clear()
      { EventChangeTracker_t::clear(); SetInputLabel(art::InputTag()); }
    
    /// Set a new event and data product label as current
    void set(art::Event const& evt, art::InputTag const& label)
      { EventChangeTracker_t::set(evt); SetInputLabel(label); }
    
    /// Update to a new data product, return true if it has changed
    bool update(DataProductChangeTracker_t const& new_prod)
      {
        if (same(new_prod)) return false;
        *this = new_prod;
        return true;
      }
    
    /// @}
    
    /// Returns a string representation of event and data product label
    operator std::string() const
      {
        return EventChangeTracker_t::operator std::string()
          + " I{" + inputLabel().encode() + "}";
      }
    
      private:
    struct LocalState_t {
      art::InputTag input_label;
    }; // LocalState_t
    
    LocalState_t state;
    
    void SetInputLabel(art::InputTag const& label)
      { state.input_label = label; }
    
  }; // DataProductChangeTracker_t
  
  
  template <typename Stream>
  decltype(auto) operator<<
    (Stream&& out, DataProductChangeTracker_t const& trk)
    { out << std::string(trk); return std::forward<Stream>(out); }
  
  
  
  /** **************************************************************************
   * @brief Detects the presence of a new event, data product or wire plane
   *
   * The state of this class describes the current data product input as the
   * event it belongs to (by its ID) and the input tag, and the current wire
   * plane in the TPC.
   */
  class PlaneDataChangeTracker_t: private DataProductChangeTracker_t {
      public:
    
    /// Default constructor: no current plane data
    PlaneDataChangeTracker_t() = default;
    
    /// Constructor: specifies current data product and TPC plane
    PlaneDataChangeTracker_t(
      art::Event const& evt, art::InputTag const& label, geo::PlaneID const& pid
      ):
      DataProductChangeTracker_t(evt, label), state{pid}
      {}
    
    /// Constructor: specifies current data product and TPC plane
    PlaneDataChangeTracker_t(
      art::EventID const& evt_id, art::InputTag const& label,
      geo::PlaneID const& pid
      ):
      DataProductChangeTracker_t(evt_id, label), state{pid}
      {}
    
    
    /// @name State query
    /// @{
    /// Returns the current input label
    using DataProductChangeTracker_t::inputLabel;
    
    /// Returns the current plane ID
    geo::PlaneID const& planeID() const { return state.plane_id; }
    
    /// Returns whether we are in the same event (the rest could differ)
    bool sameProduct(PlaneDataChangeTracker_t const& as) const
      { return DataProductChangeTracker_t::same(as); }
    
    /// Returns whether we have the same data product and TPC as "as"
    bool sameTPC(PlaneDataChangeTracker_t const& as) const
      {
        return sameEvent(as)
          && (static_cast<geo::TPCID const&>(planeID()) == as.planeID());
      }
    
    /// Returns whether we have the same plane data as "as"
    bool same(PlaneDataChangeTracker_t const& as) const
      { return sameEvent(as) && (planeID() == as.planeID()); }
    
    /// Returns whether there is a data product and plane
    bool isValid() const
      { return DataProductChangeTracker_t::isValid() && planeID().isValid; }
    
    /// Returns whether data product and TPC plane are the same as in "as"
    bool operator== (PlaneDataChangeTracker_t const& as) const
      { return same(as); }
    
    /// Returns whether data product or TPC plane are different than in "than"
    bool operator!= (PlaneDataChangeTracker_t const& than) const
      { return !same(than); }
    /// @}
    
    
    /// @name State change
    /// @{
    /// Forget the current data product
    void clear()
      { DataProductChangeTracker_t::clear(); SetPlaneID(geo::PlaneID()); }
    
    /// Set a new event and data product label as current
    void set(
      art::Event const& evt, art::InputTag const& label, geo::PlaneID const& pid
      )
      { DataProductChangeTracker_t::set(evt, label); SetPlaneID(pid); }
    
    /// Update to a new data product, return true if it has changed
    bool update(PlaneDataChangeTracker_t const& new_prod)
      {
        if (same(new_prod)) return false;
        *this = new_prod;
        return true;
      }
    
    /// @}
    
    /// Returns a string representation of event and data product label
    operator std::string() const
      {
        return DataProductChangeTracker_t::operator std::string()
          + " " + std::string(planeID());
      }
    
      private:
    struct LocalState_t {
      geo::PlaneID plane_id;
    }; // LocalState_t
    
    LocalState_t state;
    
    void SetPlaneID(geo::PlaneID const& pid) { state.plane_id = pid; }
    
  }; // PlaneDataChangeTracker_t
  
  
  template <typename Stream>
  decltype(auto) operator<< (Stream&& out, PlaneDataChangeTracker_t const& trk)
    { out << std::string(trk); return std::forward<Stream>(out); }
  
  
} // namespace util

#endif // UTIL_CHANGETRACKERS_H
