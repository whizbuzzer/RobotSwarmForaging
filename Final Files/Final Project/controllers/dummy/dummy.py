"""dummy controller."""

# You may need to import some classes of the controller module. Ex:
#  from controller import Robot, Motor, DistanceSensor
from controller import Robot, Emitter, Receiver

TIME_STEP = 64
# create the Robot instance.
robot = Robot()
#print(robot.getName())

# get the time step of the current world.
# timestep = int(robot.getBasicTimeStep())

# You should insert a getDevice-like function in order to get the
# instance of a device of the robot. Something like:
#  motor = robot.getMotor('motorname')
#  ds = robot.getDistanceSensor('dsname')
#  ds.enable(timestep)
receiver = robot.getReceiver("receiver")
receiver.enable(TIME_STEP)
emitter = robot.getEmitter("emitter")
emitter.setRange(0.17)
# print("dummy's emission range is:",emitter.getRange())

#print(receiver.getSamplingPeriod())

#print(receiver.getChannel())
#print(receiver.getSignalStrength())
# Main loop:
# - perform simulation steps until Webots is stopping the controller
while robot.step(TIME_STEP) != -1:
#    print(receiver.getQueueLength())
    if receiver.getQueueLength() > 0:
        message = receiver.getData()
        print(message)
        # print(receiver.getQueueLength())
        receiver.nextPacket()
        
    message = str.encode("aye yo")
    emitter.send(message)
    # Read the sensors:
    # Enter here functions to read sensor data, like:
    #  val = ds.getValue()
    #if receiver.getQueueLength() > 0:
        #message = receiver.getData()
        #print(receiver.getQueueLength())
        #receiver.nextPacket()
    # dataList = struct.unpack("chd",message)
        #if message == "food":
        #    food += 1
    # Process sensor data here.

    # Enter here functions to send actuator commands, like:
    #  motor.setPosition(10.0)
    pass

# Enter here exit cleanup code.
