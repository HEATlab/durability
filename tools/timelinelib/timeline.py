##
# \file timeline.py
# \brief Timeline data structure for all simulations.
# \details This data structure replaces the old legacy simulation code
#   used by robot brunch in the past. It allows the conversion of an
#   STN to a linear timeline, and also allows converting this timeline
#   back to an STN.
#

import os, sys
import random

# use CURDIR to find tools directory, which we then add to the python path.
CURDIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(CURDIR+'/../')

import distgenlib
from stntools import STN, Vertex, Edge, Triangle

## Do the debuggy thingies
DEBUG = False
UNIT_CONVERSION = 1000

class TimelineEvent:
  ##
  # \fn __init__
  #
  # @param vertex STN Vertex object to use for the Timepoint.
  # @param min_time The minimum time this TimelineEvent could occur.
  # @param max_time The maximum time this TimelineEvent could occur.
  # @param contingent_parent TimelineEvent representing a contingent edge
  #   from the parent to this TimelineEvent(self). Leave at None for
  #   requirement edges.
  # @param children List of TimelineEvents that depend on this TimelineEvent.
  # @param parents List of TimelineEvents that this TimelineEvent depends on.
  def __init__(self,
               vertex,
               min_time,
               max_time,
               children,
               parents,
               contingent_parent):

    self.vertex = vertex
    self.min_time = min_time
    self.max_time = max_time
    self.children = children
    self.parents = parents

    self.executed = False # Has this TimelineEvent been executed yet?
    self.time = None # The actual time that this executed event occured at.
    self.placed = False # Has this time point been placed already?

    # Does this TimelineEvent have a parent that is connected via a contingent
    # edge? This is what the contingent_parent variable represents.
    self.contingent_parent = contingent_parent
    # If this is a contingent edge, then we'll need to store the sampled
    # time it's offset from its parent.
    self.contingent_offset = None
    # If this is a deterministic edge, we want to store the minimum time 
    # for that edge.  This is in place so that events do not get executed
    # before it is possible for them to be executed.
    self.deterministic_offset = []
    # String representing which .samples file to use for the contingent edge.
    self.contingent_dist = None


  ## String representation
  def __str__(self):
    return ("TimelineEvent(vertex={0}, min_time={1}, max_time={2}"
            +", children={3}, parents={4}, contingent_parent{5})"
           ).format(str(self.vertex),
                    self.min_time,
                    self.max_time,
                    self.children,
                    self.parents,
                    self.contingent_parent)

  ##
  # \fn set_child(child)
  # \brief adds child to the children list of parent, and parent to the
  #  parents list of child.
  #
  # @post children list of parent is updated with child, and parents list
  #  of child is updated with parent. Note, this modifies the child as well.
  def set_child(self, child):
    self.children.append(child)
    child.parents.append(self)

# =============================================================================
# T I M E L I N E
# =============================================================================

class Timeline:
  def __init__(self):

    ## Unsorted list of events that still need to be placed.
    # NOTE: May be removed?
    self.unplaced_events = []
    ## The events that have been placed on the timeline, sorted by when
    # they will occur (earliest first)
    self.placed_events = []
    ## Current time in the timeline.
    self.current_time = 0.0
    ## Number of events the timeline has passed
    self.events_passed = 0
    ## The internal samples map used when resampling contingent edges.
    self.samp_map = None


  ##
  # \fn advance
  # \brief Advances the timeline exactly one event.
  # \details Finds the next available timepoint event, and executes it.
  #   Then updates the timeline if any of the future children can be
  #   placed.
  #
  # @return Returns either a tuple of the form (Vertex, time, contingent_bool)
  #   of the executed timepoint, or returns None if the timeline is empty.
  # @post The timeline advances one event (timepoint), if the next event
  #   exists.
  def advance(self):

    # If we've passed all the placed events already, we've reached the end
    # of the timeline. Stop here, and return None
    if self.events_passed >= len(self.placed_events):
      return None

    # We still have more events left to iterate through, use these!
    next_event = self.placed_events[self.events_passed]
    next_event.executed = True
    # May seem odd, but min_time could have changed since the last time
    # we set it. This way, we can nail down the actual time the event occured.
    next_event.time = next_event.min_time
    # Move the timeline forwards
    self.current_time = next_event.time
    # Increment the number of events we've passed in the timeline.
    self.events_passed += 1

    # Check all the dependents of our next event,
    # see if any of them can be placed on the timeline.
    for child in next_event.children:

      # By default, assume the child can be placed.
      # Then later on contradict that statement.
      child_can_be_placed = True

      if child.contingent_parent == None:
        # The child only has requirement edges connecting to it.

        # -- Begin inner for loop------
        for parent in child.parents:
          if not parent.executed:
            child_can_be_placed = False
            break
        # -- End inner for loop--------

      else:
        # The child has a contingent edge connecting to it somehow.
        if not child.contingent_parent.executed:
          child_can_be_placed = False

      if child.placed:
        child_can_be_placed = False

      if child_can_be_placed:
        # All the dependencies of the child have been executed
        # We can thus place its minimum time on the timeline.
        self.place_event(child)
        # Indicates if one of the children was placed

    # --------End outer for loop----------

    tup = (next_event.vertex,
           next_event.time,
           (next_event.contingent_parent != None))
    return tup


  ##
  # \fn peek
  # \brief View the first item on the timeline, without executing.
  # \details Very similar to peek() in other data structures,
  #   this allows the user to view the front of the timeline without
  #   changing anything.
  # @return Returns either a tuple of the form (Vertex, time, contingent_bool) of
  #   the executed timepoint, or returns None if the timeline is empty.
  def peek(self):
    next_event = self.placed_events[self.events_passed]
    tup = (next_event.vertex,
           next_event.min_time,
           (next_event.contingent_parent != None))
    return tup


  ##
  # \fn place_event
  # \brief Helper function for advance() to allow implementation of insert
  #   sort.
  # \details Inserts a TimelineEvent into the timeline in the correct order.
  #   In other words, it places an unplaced event in the right order.
  # \note May modify the event you insert, as the min time will be forced to
  #   be at least the current time. Thus, using the event after inserting
  #   is now invalid for further use.
  #
  # @param event Event to insert into the timeline.
  # @return a boolean value, False if min_time > max_time
  def place_event(self, event):
    # Check if we're dealing with a contingent event.
    if event.contingent_parent == None:
      earliest_time = self.current_time
      for tup in event.deterministic_offset:
        time = tup[0].time + tup[1]
        if time > earliest_time:
          earliest_time = time

      # We can't do an event that's earlier than now!
      if event.min_time < earliest_time:
        # Modify the event's min time and set it to now (!WARNING!)
        event.min_time = earliest_time

      if (event.min_time > event.max_time and DEBUG):
        print(("Node {0}'s min time ({1}) was greater than its max time"
          +" ({2})").format(event.vertex.nodeID,
                            event.min_time,
                            event.max_time))
    else:
      # We need to handle contingent edges differently, as we know exactly
      # when they should occur.

      if not event.contingent_parent.executed:
        # Events with a contingent edge leading to them must have their
        # single parent executed to be able to be placed.
        raise RuntimeError(("Attempted to place contingent event {0}"
            +" without executing parent {1}.").format(str(e),
            str(contingent_parent)))
      else:
        # We can't set the actual time, as we haven't executed the contingent
        # edge. So instead set the min and max time to the same value.
        event.min_time = event.contingent_parent.time+event.contingent_offset

        if event.min_time > event.max_time:
          consistent = False

        event.max_time = event.min_time

    # Loop flag
    successfully_inserted = False
    for i in range(self.events_passed,len(self.placed_events)):
      # Iterate through the remaining events until we find one that is after
      # our newly inserted event.
      if event.min_time < self.placed_events[i].min_time:
        self.placed_events.insert(i,event)
        event.placed = True
        successfully_inserted = True
        break
    if not successfully_inserted:
      # We couldn't insert because nothing in the list was before ours,
      # so just append instead.
      self.placed_events.append(event)
      event.placed = True


  ##
  # \fn loadtimeline(self, stn)
  # \brief takes in an STN object and creates Timeline events from each
  #  timepoint in the STN.  All of these events are put in
  #  unplacedEvents, and the zero timepoint is put in placedEvents.
  # @param stn an STN object (not a json)
  # @return nothing
  # @post unplacedEvents is filled with all TimelineEvents except for
  #  the zero timepoint which is in placedEvents.  The zero timepoint
  #  event is given a min_time of 0 and a max_time of infinity.
  def loadtimeline(self, stn):
    # parents and childrens dictionaries are meant to hold info about
    # parent and childrens until all timeline events have been updated.
    event_dict = {}
    new_parents = {}
    new_children = {}

    # Initializes a dictionary of form {vertex:TimelineEvent} to store
    # timeline events while they are being modified.
    for vertexID in stn.verts:
      event_dict[vertexID] = TimelineEvent(vertex=stn.verts[vertexID],
                                         min_time=None,
                                         max_time=None,
                                         children=[],
                                         parents=[],
                                         contingent_parent=None)
      new_parents[vertexID] = []
      new_children[vertexID] = []

    # Initializes special properties of zero timepoint
    zero_id = 0
    event_dict[zero_id].min_time = 0.0
    event_dict[zero_id].max_time = "inf"
    event_dict[zero_id].contingent_parent = None

    for edge in stn.edges:

      edgeOb = stn.edges[edge] # The actual edge object from the STN.
      vert1 = stn.edges[edge].i # What the edge is pointing *from*.
      vert2 = stn.edges[edge].j # What the edge is pointing *to*.

      # Checks to see if it is a zero edge or not. If it is a zero edge
      # it can be used to update the min and max time.  If it is a regular
      # edge it can be used to determine children/parents and incoming edge
      # type.
      if vert1 == zero_id:

        event_dict[vert2].min_time = -edgeOb.Cji
        if edgeOb.Cij == None:
          event_dict[vert2].max_time = "inf"
        else:
          event_dict[vert2].max_time = edgeOb.Cij

      # Otherwise it is an edge between two regular vertices so it
      # determines if it is a contingent edge or not.
      else:

        if edgeOb.distribution != None:
          # Set the parent to a TimelineEvent.
          parent_event = event_dict[edgeOb.i]
          event_dict[vert2].contingent_parent = parent_event
          event_dict[vert2].contingent_dist = edgeOb.distribution
        else:
          # Stores a tuple of (parent, minimum edge time)
          event_dict[vert2].deterministic_offset += [(event_dict[edgeOb.i], -edgeOb.Cji)]

        # Updates the children and parents. Doesn't do this if the edge
        # leads from the zero timepoint otherwise every timeline event
        # would be set as a child of the zero node.
        new_parents[vert2] += [vert1]
        new_children[vert1] += [vert2]

    # Places all TimelineEvents in their proper lists
    for key,value in event_dict.iteritems():
      # Initializes the zero node timeline events with the proper
      # children/parents.  This ensures that place_events starts out
      # placing the events that occur first chronologically.
      if new_parents[key] == [] and key != 0:
        new_parents[key] += [0]
        event_dict[0].children += [event_dict[key]]

      # Updates timeline events with the proper children and parents.
      for vert in new_parents[key]:
        value.parents += [event_dict[vert]]

      for vert in new_children[key]:
        value.children += [event_dict[vert]]

      # Puts the zero timepoint in placed_events and all other timeline
      # events in unplaced_events.
      if key == 0:
        self.place_event(value)
      else:
        self.unplaced_events += [value]


  ##
  # \fn update_timeline(stn)
  # \brief uses an stn to update all the events on the timeline
  # \details updates the min and max time of each timeline event to reflect
  #  any changes that an STN may have undergone.  This is intended for use
  #  in DREA but not SREA or early first. This is to be used in conjunction
  #  with update_stn to communicate between the timeline and an STN.
  #
  # @param stn an STN object
  # @post min and max time of timeline events in unplaced_events are updated.
  def update_timeline(self, stn):
    for event in self.unplaced_events:
      # We don't want to update timepoints that have incoming contingent edges
      # since the time they occur will be decided by the execution time of the
      # preceding timepoint.
      if event.contingent_offset == None and event.executed == False:
        vert = event.vertex
        edge = stn.edges[(0, vert.nodeID)]
        
        new_min_time = -edge.Cji
        new_max_time = edge.Cij
        # Makes sure that no time travelling occurs
        if new_min_time < self.current_time:
          new_min_time = self.current_time
        if new_max_time < self.current_time:
          new_max_time = self.current_time
        # Gets min and max times from stn and transfers them to timeline.
        event.min_time = new_min_time
        event.max_time = new_max_time


  ##
  # \fn update_stn(stn)
  # \brief updates the zero edges of stn with executed timeline events.
  # \details loops through all the executed events and updates their corresponding
  #  vertices in the stn. This is to be used in conjunction with update_timeline
  #  to communicate between the timeline and the STN.
  #
  # @param stn an STN object
  # @return an updated STN object
  # @post corresponding min and max times in the STN are modified.
  def updated_stn(self, stn):
    copy = stn.copy()

    for index in range(len(self.placed_events)):
      if index != 0:
        event = self.placed_events[index]
        vert = event.vertex
        if event.time != None:
          copy.edges[(0, vert.nodeID)].Cij = int(event.time)
          copy.edges[(0, vert.nodeID)].Cji = -int(event.time)

    return copy

  ##
  # \fn check_stn(stn)
  # \brief checks the timeline against an stn
  # \details takes the given stn and makes sure that all executed timeline
  #  events occur at valid times in the stn.
  #
  # @param stn an STN object
  # @return a boolean value indicating if the check was successful
  def check_stn(self, stn):
    # Boolean indicating if the timeline is valid for the given stn.
    success = True
    # Checks all executed events to see if they fit in the stn.
    for event in self.placed_events:

      vert = event.vertex
      if event.executed == True and vert.nodeID != 0:

        zero_edge = stn.edges[(0,vert.nodeID)]
        upper_limit = zero_edge.Cij
        lower_limit = -zero_edge.Cji
        # Checks if the execution time for the timeline event is outside the
        # bounds of the corresponding vertex
        if event.time < lower_limit or event.time > upper_limit:
          success = False

    return success



  ##
  # \fn set_sample_map
  # \brief Set the sample map to sample for calculating contingent edge
  #   durations
  #
  # @param samp_map Sample map to use when resampling contingent edges.
  # @post Any further resample_contingent calls will use the new sample map set
  #   here.
  def set_sample_map(self, samp_map):
    self.samp_map = samp_map


  ##
  # \fn resample_contigent
  # \brief Resamples the contingent edges in the loaded STN.
  # \details Using the set invcdf map (see set_invcdfmap), resample the
  #   contingent edges so that they are ready for a different simulation.
  #
  # @post Changes all the internal contingent edge values to new randomly
  #   sampled values.
  def resample_contingent(self):
    for e in self.unplaced_events:
      if e.contingent_parent != None:
        e.contingent_offset = random.choice(
            self.samp_map[e.contingent_dist]) * UNIT_CONVERSION


  ##
  # \fn get_executed()
  # \brief displays timelineEvents that have been executed, meaning that
  #  a specific time has been picked for those events.
  # @return a dictionary of the form {vertex:time}
  def get_executed(self):
    executed = {}
    # This only loops through placedEvents because events can only be executed
    # if they are placed.
    for event in self.placedEvents:
      if event.executed == True:
        executed[event.vertex] = event.time
    return executed


  ##
  # \fn get_timeline()
  # \brief displays timelineEvents that have been placed on the timeline
  # @return a dictionary of the form {vertex:time}
  def get_timeline(self):
    timeline = {}
    for event in self.placedEvents:
      timeline[event.vertex] = event.time
    return timeline


  ##
  # \fn zedge_executed
  # \brief Returns Z-Edge represenation of all the executed timepoints.
  #
  # @return Returns a List of Edges C0i which determine when timepoints
  #   occur relative to the z-timepoint (requirement edges with a single
  #   viable time).
  def zedge_executed(self):
    edge_list = []
    for event in self.placedEvents:
      if event.executed == True:
        edge_list.append(Edge(0, event.vertex, event.time, event.time))
    return edge_list


