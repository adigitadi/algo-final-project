import { ChangeEventHandler, useState } from "react";
import { Rider, Taxi } from "../../types";

type Props = {
  height: number;
  width: number;
  riders: Rider[];
  numRiders: number;
  taxis: Taxi[];
  numTaxis: number;
  onHeightChange: ChangeEventHandler<HTMLInputElement> | undefined;
  onWidthChange: ChangeEventHandler<HTMLInputElement> | undefined;
  onNumRidersChange: ChangeEventHandler<HTMLInputElement> | undefined;
  onNumTaxisChange: ChangeEventHandler<HTMLInputElement> | undefined;
  populateRiders: (event: React.MouseEvent<HTMLButtonElement>) => void;
  clearRiders: (event: React.MouseEvent<HTMLButtonElement>) => void;
  populateTaxis: (event: React.MouseEvent<HTMLButtonElement>) => void;
  clearTaxis: (event: React.MouseEvent<HTMLButtonElement>) => void;
  calculatePaths: (event: React.MouseEvent<HTMLButtonElement>) => void;
};

const RiderList = (props: Props) => {
  const [showNewRiderCard, setShowNewRiderCard] = useState(false);

  return (
    <div className="h-screen max-h-screen bg-slate-800 z-50 p-2  overflow-y-scroll">
      <div className=" justify-center pt-2 pb-2 ">
        <div className="rounded-lg shadow-lg bg-white max-w-sm ">
          <div className="p-3">
            <div className="flex justify-center">
              <div className="mb-0 xl:w-96">
                <button
                  type="button"
                  onClick={props.calculatePaths}
                  className="inline-block px-6 py-2.5 h-20 w-full bg-green-500 text-white font-medium text-xs leading-tight uppercase rounded shadow-md hover:bg-green-600 hover:shadow-lg focus:bg-green-600 focus:shadow-lg focus:outline-none focus:ring-0 active:bg-green-700 active:shadow-lg transition duration-150 ease-in-out"
                >
                  Calculate Paths
                </button>
              </div>
            </div>
          </div>
        </div>
      </div>

      <div className=" justify-center ">
        <div className="rounded-lg shadow-lg bg-white max-w-sm ">
          <div className="p-3">
            <h5 className="text-gray-900 text-xl font-medium mb-2 ">
              {props.riders.length ? (
                <>Riders waiting...</>
              ) : (
                <>No Riders Yet.</>
              )}
            </h5>
            <>
              <div className="h-36 overflow-y-auto b-black">
                {props.riders.map((rider) => {
                  return (
                    <div className="pb-5" key={rider.id}>
                      <h5 className="text-left text-gray-900 text-md font-small ">
                        {rider.name}
                      </h5>
                      <p className="text-left text-gray-900 text-md font-small ">
                        Position: Node {rider.source}
                      </p>
                    </div>
                  );
                })}
              </div>
            </>
            <div className="grid grid grid-cols-1 gap-4">
              {/* TODO. Only let users populate list once. Geo codes will stack. */}
              <div className="flex justify-center">
                <div className="mb-0 xl:w-96">
                  <label className="form-label inline-block mb-2 text-gray-700">
                    <b>Add # Riders</b>
                  </label>
                  <input
                    type="number"
                    className="
                      form-control
                      block
                      w-full
                      px-3
                      py-1.5
                      text-base
                      font-normal
                      text-gray-700
                      bg-white bg-clip-padding
                      border border-solid border-gray-300
                      rounded
                      transition
                      ease-in-out
                      m-0
                      focus:text-gray-700 focus:bg-white focus:border-blue-600 focus:outline-none
                    "
                    id="riderInput"
                    placeholder="0"
                    onChange={props.onNumRidersChange}
                    value={props.numRiders}
                  />
                </div>
              </div>
              <button
                onClick={props.populateRiders}
                type="button"
                className=" inline-block px-6 py-2.5 bg-blue-600 text-white font-medium text-xs leading-tight uppercase rounded shadow-md hover:bg-blue-700 hover:shadow-lg focus:bg-blue-700 focus:shadow-lg focus:outline-none focus:ring-0 active:bg-blue-800 active:shadow-lg transition duration-150 ease-in-out"
              >
                Populate Riders
              </button>
              <button
                onClick={props.clearRiders}
                type="button"
                className=" inline-block px-6 py-2.5 bg-red-600 text-white font-medium text-xs leading-tight uppercase rounded shadow-md hover:bg-red-700 hover:shadow-lg focus:bg-red-700 focus:shadow-lg focus:outline-none focus:ring-0 active:bg-red-800 active:shadow-lg transition duration-150 ease-in-out"
              >
                Clear Riders
              </button>
            </div>
          </div>
        </div>
      </div>
      <div className=" justify-center pt-2">
        <div className="rounded-lg shadow-lg bg-white max-w-sm ">
          <div className="p-3">
            <h5 className="text-gray-900 text-xl font-medium mb-2 ">
              {props.taxis.length ? <>Taxis waiting...</> : <>No Taxis Yet.</>}
            </h5>
            <>
              <div className="h-36 overflow-y-auto b-black">
                {props.taxis.map((taxi) => {
                  return (
                    <div className="pb-5" key={taxi.id}>
                      <h5 className="text-left text-gray-900 text-md font-small ">
                        {taxi.name}
                      </h5>
                      <p className="text-left text-gray-900 text-md font-small ">
                        Position: Node {taxi.source}
                      </p>
                      <p className="text-left text-gray-900 text-md font-small ">
                        Capacity: {taxi.capacity}
                      </p>
                    </div>
                  );
                })}
              </div>
            </>

            <div className="grid grid grid-cols-1 gap-4">
              <div className="flex justify-center">
                <div className="mb-0 xl:w-96">
                  <label className="form-label inline-block mb-2 text-gray-700 ">
                    <b>Add # Taxis</b>
                  </label>
                  <input
                    type="number"
                    className="
                      form-control
                      block
                      w-full
                      px-3
                      py-1.5
                      text-base
                      font-normal
                      text-gray-700
                      bg-white bg-clip-padding
                      border border-solid border-gray-300
                      rounded
                      transition
                      ease-in-out
                      m-0
                      focus:text-gray-700 focus:bg-white focus:border-blue-600 focus:outline-none
                    "
                    id="riderInput"
                    placeholder="0"
                    onChange={props.onNumTaxisChange}
                    value={props.numTaxis}
                  />
                </div>
              </div>

              {/* TODO. Only let users populate list once. Geo codes will stack. */}
              <div className="flex justify-center"></div>
              <button
                onClick={props.populateTaxis}
                type="button"
                className=" inline-block px-6 py-2.5 bg-yellow-600 text-white font-medium text-xs leading-tight uppercase rounded shadow-md hover:bg-yellow-700 hover:shadow-lg focus:bg-yellow-700 focus:shadow-lg focus:outline-none focus:ring-0 active:bg-yellow-800 active:shadow-lg transition duration-150 ease-in-out"
              >
                Populate Taxis
              </button>
              <button
                onClick={props.clearTaxis}
                type="button"
                className=" inline-block px-6 py-2.5 bg-red-600 text-white font-medium text-xs leading-tight uppercase rounded shadow-md hover:bg-red-700 hover:shadow-lg focus:bg-red-700 focus:shadow-lg focus:outline-none focus:ring-0 active:bg-red-800 active:shadow-lg transition duration-150 ease-in-out"
              >
                Clear Taxis
              </button>
            </div>
          </div>
        </div>
      </div>
      <div className=" justify-center pt-2 pb-2 ">
        <div className="rounded-lg shadow-lg bg-white max-w-sm ">
          <div className="p-3">
            <h5 className="text-gray-900 text-xl font-medium mb-2 ">Graph</h5>
            <div className="flex justify-center">
              <div className="mb-0 xl:w-96">
                <label className="form-label inline-block mb-2 text-gray-700">
                  Height
                </label>
                <input
                  type="number"
                  className="
                      form-control
                      block
                      w-full
                      px-3
                      py-1.5
                      text-base
                      font-normal
                      text-gray-700
                      bg-white bg-clip-padding
                      border border-solid border-gray-300
                      rounded
                      transition
                      ease-in-out
                      m-0
                      focus:text-gray-700 focus:bg-white focus:border-blue-600 focus:outline-none
                    "
                  id="riderInput"
                  placeholder="0"
                  onChange={props.onHeightChange}
                  value={props.height as number}
                />
              </div>
            </div>
            {/* <div className="flex justify-center">
              <div className="mb-0 xl:w-96">
                <label className="form-label inline-block mb-2 text-gray-700">
                  Width
                </label>
                <input
                  type="number"
                  className="
                      form-control
                      block
                      w-full
                      px-3
                      py-1.5
                      text-base
                      font-normal
                      text-gray-700
                      bg-white bg-clip-padding
                      border border-solid border-gray-300
                      rounded
                      transition
                      ease-in-out
                      m-0
                      focus:text-gray-700 focus:bg-white focus:border-blue-600 focus:outline-none
                    "
                  id="riderInput"
                  placeholder="0"
                  onChange={props.onWidthChange}
                  value={props.width as number}
                />
              </div>
            </div> */}
          </div>
        </div>
      </div>
    </div>
  );
};

export default RiderList;
