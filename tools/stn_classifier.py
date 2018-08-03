import json
import os
from sklearn import svm
import stnConverter
import randomSchedule
import torch
import torch.nn.functional as F

with open(os.path.abspath('./Data/' + 'ALL_OF_THE_DATA_30000' + '.json'), 'r') as fp:
	datastuff = json.load(fp)

twoData = []
newData = []
walks = []

for i in range(30000):
	if datastuff[1][i] == 5:
		twoData.append(stnConverter.dictToMatrix(randomSchedule.fixStupidJSON(datastuff[0][i])))
		walks.append(datastuff[6][i])

for mat in twoData:
	newData.append(mat[0]+mat[1]+mat[2]+mat[3]+mat[4]+mat[5])

lengthTrain = int(0.8 * len(newData))

trainingSet = newData[:lengthTrain]
trainingLabels = walks[:lengthTrain]

testSet = newData[lengthTrain:]
testLabels = walks[lengthTrain:]

class Net(torch.nn.Module):
    def __init__(self, n_feature, n_hidden, n_output):
        super(Net, self).__init__()
        self.hidden = torch.nn.Linear(n_feature, n_hidden)   # hidden layer
        # self.hidden1 = torch.nn.Linear(n_hidden, n_hidden)   # hidden layer
        # self.hidden2 = torch.nn.Linear(n_hidden, n_hidden)   # hidden layer
        # self.hidden3 = torch.nn.Linear(n_hidden, n_hidden)   # hidden layer
        # self.hidden4 = torch.nn.Linear(n_hidden, n_hidden)   # hidden layer
        # self.hidden5 = torch.nn.Linear(n_hidden, n_hidden)   # hidden layer
        # self.hidden6 = torch.nn.Linear(n_hidden, n_hidden)   # hidden layer
        # self.hidden7 = torch.nn.Linear(n_hidden, n_hidden)   # hidden layer
        # self.hidden8 = torch.nn.Linear(n_hidden, n_hidden)   # hidden layer
        self.out = torch.nn.Linear(n_hidden, n_output)   # output layer

    def forward(self, x):
        x = F.relu(self.hidden(x))      # activation function for hidden layer
        # x = F.relu(self.hidden1(x))
        # x = F.relu(self.hidden2(x))
        # x = F.relu(self.hidden3(x))
        # x = F.relu(self.hidden4(x))
        # x = F.relu(self.hidden5(x))
        # x = F.relu(self.hidden6(x))
        # x = F.relu(self.hidden7(x))
        # x = F.relu(self.hidden8(x))
        x = self.out(x)
        return x

net = Net(n_feature=36, n_hidden=50, n_output=1)

inp = torch.tensor(trainingSet)
out = torch.tensor(trainingLabels)
out = out.unsqueeze(1)



optimizer = torch.optim.SGD(net.parameters(), lr=0.02)
# loss_func = torch.nn.ReLU()  # the target label is NOT an one-hotted
loss_func = torch.nn.L1Loss()

for t in range(10000):
    # Forward pass: Compute predicted y by passing x to the model
    y_pred = net(inp)
    loss = loss_func(y_pred, out)
    print(t, loss.item())

    # Zero gradients, perform a backward pass, and update the weights.
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()

PATH = os.path.abspath('./Data/' + 'CHEB_RANDOM_WALK_PREDICTINATOR' + '.json')

# net.load_state_dict(torch.load(PATH))

success = 0
testingTest = 0

if testingTest == 1:
	trainingSet = testSet
	trainingLabels = testLabels

test = torch.tensor(trainingSet)
preds = net(test)
l1 = []
for i in range(len(preds)):
	l1.append(preds[i][0])

margins = randomSchedule.listMargin(l1, trainingLabels)
for marg in margins:
	if abs(marg) < 0.2:
		success += 1

# print preds
# print testLabels
print l1
print trainingLabels
print success / (len(trainingSet) + 0.0)


# torch.save(net.state_dict(), PATH)

# twoData = []
# twoDataLabels = []

# for i in range(30000):
# 	chebWalk = datastuff[10][i]
# 	centWalk = datastuff[8][i]
# 	midWalk = datastuff[10][i]
# 	margin = randomSchedule.calcMargin(chebWalk, centWalk)

# 	if datastuff[1][i] == 2 and abs(margin) > 0.15:
# 		twoData.append(stnConverter.dictToMatrix(randomSchedule.fixStupidJSON(datastuff[0][i])))

# 		if chebWalk > centWalk:
# 			twoDataLabels.append(0)
# 		else:
# 			twoDataLabels.append(1)
# # print len(twoData)
# newData = []

# for mat in twoData:
# 	newData.append(mat[0]+mat[1]+mat[2])

# lengthTrain = int(0.8 * len(newData))

# trainingSet = newData[:lengthTrain]
# trainingLabels =twoDataLabels[:lengthTrain]

# testSet = newData[lengthTrain:]
# testLabels = twoDataLabels[lengthTrain:]

# clf = svm.SVC()
# clf.fit(trainingSet, trainingLabels)

# success = 0

# labels = clf.predict(testSet)

# for i in range(len(labels)):
# 	if labels[i] == testLabels[i]:
# 		success += 1
# print labels
# print testLabels
# print success / (len(testSet) + 0.0)
# print clf.kernel