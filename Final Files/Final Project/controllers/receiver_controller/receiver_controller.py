"""receiver_controller controller."""

# You may need to import some classes of the controller module. Ex:
#  from controller import Robot, Motor, DistanceSensor
from controller import Robot, Emitter, Receiver

# get the time step of the current world.
TIME_STEP = 64

# create the Robot instance.
robot = Robot()

emitter = robot.getEmitter("nest_emitter")

# food_received = 0

receiver = robot.getReceiver('nest_receiver')
receiver.enable(TIME_STEP)

food_received = 0  # food flag
food_count = 0  # counting amount of food returned

# Main loop:
# - perform simulation steps until Webots is stopping the controller
while robot.step(TIME_STEP) != -1:
 
#    print(receiver.getQueueLength())

    place_food_message = str.encode("nest")
    emitter.send(place_food_message)
    
    if receiver.getQueueLength() > 0:
        message_received = receiver.getData()
    
        message_received = message_received.decode("utf-8")
        #print(message_received)
      
        #print("message_received  ]==  "+message_received)    
        receiver.nextPacket()
        
        if message_received == "food received":
            food_received += 1
            print("no.of food:",food_received)   
#    pass

# Enter here exit cleanup code.
